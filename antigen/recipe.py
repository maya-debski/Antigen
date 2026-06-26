from collections import defaultdict

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, List, Optional, Callable, Union
import yaml
import re
import logging

import numpy as np

logger = logging.getLogger('antigen.recipe')

@dataclass
class InputField:
    """Represents a single input field with its validation rules."""
    type: str
    required: bool = False
    description: str = ""
    units: str = ""
    default: Any = None
    validate: Dict[str, Any] = field(default_factory=dict)
    source: Optional[str] = None

    def validate_value(self, value: Any) -> Optional[str]:
        """
        Validates a given value against various conditions such as type, required
        rules, and additional validation constraints. This method ensures that the
        provided value meets specified requirements or raises appropriate validation
        messages.

        Args:
            value: The value to be validated. Can be of any type depending on the
                validation rules defined.

        Returns:
            Optional[str]: Returns an error message string if validation fails;
                otherwise, returns None.
        """
        if value is None:
            return "Required field is missing" if self.required else None

        # Type validation
        if not self._validate_type(value):
            return f"Expected type {self.type}, got {type(value).__name__}"

        # Validation rules
        for rule, rule_value in self.validate.items():
            if rule == 'min_length' and len(value) < rule_value:
                return f"Length must be at least {rule_value}"
            elif rule == 'pattern' and not re.match(rule_value, str(value)):
                return f"Must match pattern {rule_value}"
            elif rule == 'choices' and value not in rule_value:
                choices = self.validate['choices']
                if isinstance(value, str):
                    if value.lower() not in [c.lower() for c in choices]:
                        return f"must be one of {choices}"
                else:
                    return f"Must be one of {rule_value}"
            elif rule == 'min' and float(value) < float(rule_value):
                return f"Must be at least {rule_value}"

        return None

    def _validate_type(self, value: Any) -> bool:
        try:
            if self.type == 'str':
                return isinstance(value, str)
            elif self.type == 'float':
                float(value)
                return True
            elif self.type == 'int':
                return isinstance(value, int)
            elif self.type == 'path':
                Path(value)
                return True
            elif self.type.startswith('list['):
                return isinstance(value, (list, tuple))
            elif self.type == 'dict':
                return isinstance(value, dict)
            elif self.type == 'ndarray':
                return isinstance(value, np.ndarray)
            return True
        except (ValueError, TypeError):
            return False

@dataclass
class FunctionExecutor:
    """Handles function loading, argument resolution, and execution."""
    function_path: str
    function_args: List[str]
    params: Dict[str, InputField]

    def build_function(self) -> Callable:
        """Import and return the function specified by function_path."""
        if not self.function_path:
            raise ValueError("No function implementation provided")

        try:
            module_path, func_name = self.function_path.rsplit('.', 1)
            module = __import__(module_path, fromlist=[func_name])
            return getattr(module, func_name)
        except (ImportError, AttributeError) as e:
            raise ImportError(f"Failed to import {self.function_path}: {str(e)}")

    def resolve_param_value(self, param_name: str, param_def: InputField, state: Dict[str, Any]) -> Any:
        """
        Resolves the value of a parameter either from a predefined default, state context,
        or specified source input.

        Args:
            param_name: The name of the parameter for which the value is being resolved.
            param_def: The input field definition which contains sourcing directions and
                defaults for resolving the parameter.
            state: A dictionary holding the current state and input mappings
                from which the parameter's value can be resolved.

        Returns:
            The resolved value of the parameter, determined from the default or a
            specified source mapping.

        Raises:
            ValueError: If the source defined in param_def is invalid or if a required
                parameter value cannot be found when expected.
        """
        if not param_def.source:
            return param_def.default

        parts = param_def.source.split('.')
        if parts[0] not in ('state', 'inputs'):
            raise ValueError(f"Invalid source '{parts[0]}' for parameter {param_name}")

        # Navigate to source
        value = state if parts[0] == 'state' else state['inputs']

        # Navigate through the path, handling missing keys gracefully
        try:
            for part in parts[1:]:
                value = value[part]
        except KeyError as e:
            # If the source path doesn't exist and parameter is not required, use default
            if not param_def.required:
                return param_def.default
            else:
                raise ValueError(f"Required parameter {param_name} from '{param_def.source}' is missing: {str(e)}")

        # Get actual parameter value
        if isinstance(value, dict):
            # Check if this dictionary is itself the parameter we want (like detector_dimensions)
            source_parts = param_def.source.split('.')
            if param_name == source_parts[-1]:
                # The dictionary itself is the parameter value we want
                return value
            # Otherwise, try to get param_name from inside the dictionary
            elif param_name in value:
                value = value[param_name]
            else:
                # If the dictionary exists but does not contain the key:
                # - For optional params, fall back to the parameter default
                # - For required params, raise a clear error
                if not param_def.required:
                    return param_def.default
                raise ValueError(f"Required parameter {param_name} from '{param_def.source}' is missing")

        if value is None and param_def.required:
            raise ValueError(f"Required parameter {param_name} from '{param_def.source}' is missing")

        return value

    def construct_args(self, state: Dict[str, Any]) -> List[Any]:
        """Build function arguments list from parameter definitions and state."""
        args = []
        for arg_name in self.function_args:
            param_def = self.params[arg_name]
            value = self.resolve_param_value(arg_name, param_def, state)
            args.append(value)
        return args

    def execute(self, state: Dict[str, Any]) -> Any:
        """Execute the function with resolved arguments."""
        func = self.build_function()
        args = self.construct_args(state)
        return func(*args)

@dataclass
class Step:
    """Represents a single step in the recipe."""
    name: str
    description: str
    params: Dict[str, InputField]
    needs: List[str]
    provides: List[str]
    label: Optional[str] = None
    executor: Optional[FunctionExecutor] = None
    save_executor: Optional[FunctionExecutor] = None  # Add this field

    @classmethod
    def from_yaml(cls, name: str, data: Dict[str, Any], step_params: Dict[str, Any] = None, label: Optional[str] = None) -> 'Step':
        """Create a Step instance from YAML data.
        
        Args:
            name: Name of the step
            data: YAML data defining the step
            step_params: Optional parameters provided in the step definition to override/fill templates
        """
        # First process any template parameters if provided
        step_params = step_params or {}
        
        # Helper function to substitute template variables in strings
        def substitute_templates(value: Any, params: Dict[str, Any]) -> Any:
            if isinstance(value, str):
                # Handle template substitution in strings
                for param_name, param_value in params.items():
                    template = "{" + param_name + "}"
                    if template in value:
                        value = value.replace(template, str(param_value))
                return value
            elif isinstance(value, list):
                return [substitute_templates(item, params) for item in value]
            elif isinstance(value, dict):
                return {k: substitute_templates(v, params) for k, v in value.items()}
            return value

        # Process parameter definitions
        params_def = data.get('params', {})
        
        # Apply template substitutions to parameter sources and other fields
        processed_params_def = substitute_templates(params_def, step_params)
        
        # Create InputField instances with processed definitions and override with step_params values
        params = {}
        for param_name, param_def in processed_params_def.items():
            field = InputField(**param_def)
            # If param is in step_params, set it as the default value
            if param_name in step_params:
                field.default = step_params[param_name]
                # Remove source if we're setting a direct value
                field.source = None
            params[param_name] = field

        # Process needs list with template substitutions
        needs = substitute_templates(data.get('needs', []), step_params)
        
        # Process provides list with template substitutions
        provides = substitute_templates(data.get('provides', []), step_params)

        # Create executor if algorithm is defined
        executor = None
        save_executor = None  # Initialize save_executor
        if 'algorithm' in data and 'function' in data['algorithm']:
            algorithm_data = substitute_templates(data['algorithm'], step_params)
            
            # Handle additional save operations if present
            if 'save' in algorithm_data:
                save_data = algorithm_data.pop('save')
                save_args = [substitute_templates(arg, step_params) for arg in save_data['args']]
                save_executor = FunctionExecutor(
                    function_path=save_data['function'],
                    function_args=save_args,
                    params=params
                )
            
            executor = FunctionExecutor(
                function_path=algorithm_data['function'],
                function_args=algorithm_data.get('args', []),
                params=params
            )

        return cls(
            name=name,
            description=data.get('description', ''),
            params=params,
            needs=needs,
            provides=provides,
            label=label,
            executor=executor,
            save_executor=save_executor  # Include save_executor in the return
        )

    def execute(self, state: Dict[str, Any]) -> None:
        """Execute the step and update state."""
        if not self.executor:
            raise ValueError(f"Step {self.name} has no executor")

        try:
            # Execute main operation
            result = self.executor.execute(state)

            results = result if isinstance(result, (tuple, list)) else [result]

            # Update state with results
            if len(results) == len(self.provides):
                for val, output_name in zip(results, self.provides):
                    state[output_name] = val
                    self.params[output_name] = self.create_input_field_from_value(val, self.name)
            else:
                expected = list(self.provides)
                expected_n = len(expected)
                try:
                    actual_n = len(results)
                except TypeError:
                    actual_n = 1
                missing_by_position = expected[actual_n:]  # what we still needed

                raise ValueError(
                    f"Unexpected result from {self.name}: expected {expected_n} outputs "
                    f"{expected}, but got {actual_n}. "
                    f"Missing (by position)={missing_by_position}."
                )

            # Execute save operation if present
            if self.save_executor:
                try:
                    self.save_executor.execute(state)
                except Exception as e:
                    logger.error(f"Failed to execute save operation for step {self.name}: {str(e)}")
                    raise

        except Exception as e:
            logger.error(f"Failed to execute step {self.name}: {str(e)}")
            raise

    def create_input_field_from_value(self, value: Any, name: str) -> InputField:
        """Create an InputField with appropriate type based on the value."""
        
        if isinstance(value, str):
            type_str = "str"
        elif isinstance(value, float):
            type_str = "float"
        elif isinstance(value, int):
            type_str = "int"
        elif isinstance(value, Path):
            type_str = "path"
        elif isinstance(value, (list, tuple)):
            type_str = "list"
        elif isinstance(value, dict):
            type_str = "dict"
        elif isinstance(value, np.ndarray):
            type_str = "ndarray"
        else:
            type_str = "any"

        return InputField(
            type=type_str,
            required=True,
            source="state",
            description=f"Output from step: {name}"
        )

@dataclass
class StepGroup:
    """Group of steps that can be iterated over."""
    name: str
    steps: List['Step']
    iterate_over: Optional[str] = None
    iterator_var: Optional[str] = None
    
    def get_iteration_source(self, state: Dict[str, Any]) -> List[Any]:
        """Get the list of items to iterate over from state."""
        if not self.iterate_over:
            return [None]  # Single iteration if no iterator specified
        
        # Parse the source path (e.g., "inputs.dataset.observation_files")
        parts = self.iterate_over.split('.')
        value = state
        for part in parts:
            value = value[part]
        return value

@dataclass
class Workflow:
    """Container for workflow steps and groups."""
    groups: List['StepGroup']

    @classmethod
    def from_yaml(cls, workflow_data: List[Dict[str, Any]], operations_data: Dict[str, Any]) -> 'Workflow':
        groups = []
        for item in workflow_data:
            # Process group
            group_steps = []
            for step_def in item['steps']:
                step_name = next(iter(step_def)) if isinstance(step_def, dict) else step_def
                step_label = step_def[step_name].get('label') if isinstance(step_def, dict) else None
                step_params = step_def[step_name].get('params', {}) if isinstance(step_def, dict) else {}
                step = Step.from_yaml(step_name, operations_data[step_name],
                                      step_params=step_params, label=step_label)
                group_steps.append(step)

            group = StepGroup(
                name=item['group'],
                steps=group_steps,
                iterate_over=item.get('iterate_over'),
                iterator_var=item.get('iterator_var')
            )
            groups.append(group)
        return cls(groups)


@dataclass
class Recipe:
    """Main recipe class that coordinates workflow execution."""
    name: str
    description: str
    workflow: Workflow
    input_schema: Dict[str, Dict[str, InputField]]
    outputs: List[str]

    @classmethod
    def load(cls, name: str, base_path: Path) -> 'Recipe':
        """Load recipe from YAML files."""
        recipe_file = base_path / "recipes" / f"{name}.yml"
        schema_file = base_path / "schema" / f"{name}_schema.yml"
        operations_file = base_path / "schema" / f"{name}_operations.yml"

        # Load recipe definition
        with open(recipe_file) as f:
            recipe_data = yaml.safe_load(f)

        # Load input schema
        with open(schema_file) as f:
            schema_data = yaml.safe_load(f)['schema']
            input_schema = {
                section: {
                    field_name: InputField(**field_def)
                    for field_name, field_def in fields.items()
                }
                for section, fields in schema_data.items()
            }

        # Load operations
        with open(operations_file) as f:
            operations_data = yaml.safe_load(f)

        # Create workflow from the workflow section if it exists
        workflow = Workflow.from_yaml(recipe_data['workflow'], operations_data['schema']['operations'])
        return cls(
            name=name,
            description=recipe_data.get('description', ''),
            workflow=workflow,
            input_schema=input_schema,
            outputs=recipe_data.get('outputs', [])
        )

    def run(self, inputs: Dict[str, Any], outdir: Path) -> Dict[str, Any]:
        """Execute the workflow and return outputs."""
        state = {
            'outdir': outdir,
            'inputs': inputs
        }

        logger.info(f"Starting recipe: {self.name}")

        for step in self.workflow.groups:
            self._execute_group(step, state)

        if errors := self.verify_recipe_outputs(state):
            raise ValueError("\n".join(errors))

        recipe_outputs = {out: state[out] for out in self.outputs}
        logger.info(f"Completed recipe: {self.name}")

        return recipe_outputs

    def _execute_group(self, group: 'StepGroup', state: Dict[str, Any]):
        """Execute a group of steps, possibly iterating over a collection."""
        logger.info(f"Starting group: {group.name}")

        iteration_items = group.get_iteration_source(state)
        for i, item in enumerate(iteration_items):
            logger.info(f"Group {group.name} iteration {i + 1}/{len(iteration_items)}")

            # Create iteration-specific state
            iter_state = state.copy()
            if group.iterator_var:
                iter_state[group.iterator_var] = item

            # Execute all steps in group
            for step in group.steps:
                self._execute_step(step, iter_state)

            # Update main state with iteration results
            state.update(iter_state)

        logger.info(f"Completed group: {group.name}")

    def _execute_step(self, step: 'Step', state: Dict[str, Any]):
        """Execute a single step with validation."""
        logger.info(f"Executing step: {step.name}")

        if errors := self.validate_state(state, step):
            raise ValueError(f"Step {step.name} state validation failed:\n" + "\n".join(errors))

        try:
            step.execute(state)
        except Exception as e:
            logger.error(f"Failed to execute step {step.name}: {str(e)}")
            raise

    def collect_inputs(self, args: Any, manifest: Dict[str, Any], config_dict: Dict[str, Any]) -> Dict[str, Any]:
        """
        Collects and organizes inputs from different sources (CLI, dataset, and
        configuration dictionary) according to the input schema.

        The method iterates over the defined input schema and extracts relevant
        fields from provided sources. It segregates inputs based on their
        source (`cli`, `dataset`, or `config`) and organizes them into a unified
        structure.

        Args:
            args: CLI arguments object or namespace containing user-provided
                inputs via the command line interface.
            manifest: A dictionary representing the dataset, containing key-value
                pairs of data fields.
            config_dict: A dictionary representing the configuration settings,
                typically provided from configuration files or other sources.

        Returns:
            Dict[str, Any]: A dictionary containing organized inputs from the CLI,
                dataset, and configuration source, structured according to the input
                schema.
        """
        inputs = {
            "cli": {},
            "dataset": {},
            "config": {}
        }

        # Map inputs according to schema
        for section, fields in self.input_schema.items():
            source = {
                "cli": args,
                "dataset": manifest,
                "config": config_dict
            }[section]

            for field_name in fields:
                if hasattr(source, field_name) if section == "cli" else field_name in source:
                    value = getattr(source, field_name) if section == "cli" else source[field_name]
                    inputs[section][field_name] = value
        return inputs

    def validate_inputs(self, inputs: Dict[str, Any]) -> List[str]:
        """
        Validates the inputs against the predefined input schema for each section and
        field. Ensures that the provided values conform to the validation rules
        defined in the input schema. Accumulates and returns a list of error messages
        if validation rules are violated.

        Args:
            inputs (Dict[str, Any]): A dictionary representing the input values where
                keys are section names and values are dictionaries containing field
                names and their respective values.

        Returns:
            List[str]: A list of error messages for all fields that failed validation.
            Each error message is in the format 'section.field_name: error'.
        """
        errors = []
        for section, fields in self.input_schema.items():
            section_inputs = inputs.get(section, {})
            for field_name, field in fields.items():
                if error := field.validate_value(section_inputs.get(field_name)):
                    errors.append(f"{section}.{field_name}: {error}")
        return errors

    def validate_state(self, state: Dict[str, Any], step: 'Step') -> List[str]:
        """
        Validates the provided state against the requirements of a specific step and
        returns a list of error messages if the validation fails.

        Args:
            state (Dict[str, Any]): Current execution state containing available inputs
                and corresponding values.
            step (Step): Step object that defines the requirements to be validated.

        Returns:
            List[str]: A list of error messages encountered during validation. If no
                errors are found, the list will be empty.
        """
        errors = []

        # Check dependencies
        missing = [need for need in step.needs if need not in state]
        if missing:
            errors.append(f"Missing required inputs: {', '.join(missing)}")

        return errors

    def verify_outputs(self, state: Dict[str, Any], step: Step) -> List[str]:
        """
        Verifies if the required outputs for a step are present in the state.

        This method checks whether all outputs required by the `provides` attribute of
        the given `step` are available in the `state`. If any required outputs are
        missing, they are added to the list of errors, and a detailed error message is
        constructed.

        Args:
            state (Dict[str, Any]): The current state that contains all the output data
                provided by previous steps.
            step (Step): The step object that specifies the outputs it is expected to
                provide.

        Returns:
            List[str]: A list of error messages indicating which expected outputs from
                the step are missing in the state.
        """
        errors = []

        missing_outputs = [out for out in step.provides if out not in state]
        if missing_outputs:
            errors.append(f"Step failed to provide outputs: {', '.join(missing_outputs)}")

        return errors

    def verify_recipe_outputs(self, state: Dict[str, Any]) -> List[str]:
        """
        Verifies that all required recipe outputs are present in the provided state.

        This function checks if any outputs defined in the instance are absent from
        the provided state. If missing outputs are found, an error message is added
        to the errors list, detailing which outputs are missing.

        Args:
            state (Dict[str, Any]): A dictionary representing the current state, where
                keys are output names and values are their corresponding data.

        Returns:
            List[str]: A list of error messages indicating which outputs, if any, are
                missing from the provided state.
        """
        errors = []

        missing_outputs = [out for out in self.outputs if out not in state]
        if missing_outputs:
            errors.append(f"Recipe missing required outputs: {', '.join(missing_outputs)}")

        return errors

    def plan(self) -> str:
        """Generate a formatted plan from the workflow."""
        lines = []
        step_num = 1

        for item in self.workflow.groups:
            lines.append(f"## Group: {item.name}")
            if item.iterate_over:
                lines.append(f"\nIterates over: {item.iterate_over}\n")

            for step in item.steps:
                lines.append(f"### {step_num}. {step.name}")
                if step.description:
                    lines.append(f"\n{step.description}\n")

                if step.params:
                    lines.append("\n**Parameters:**")
                    for param_name, param in step.params.items():
                        param_desc = f"- `{param_name}` ({param.type})"
                        if param.description:
                            param_desc += f": {param.description}"
                        if param.units:
                            param_desc += f" [{param.units}]"
                        if param.required:
                            param_desc += " *(required)*"
                        lines.append(param_desc)

                if step.provides:
                    lines.append("\n**Provides:**")
                    for provide in step.provides:
                        lines.append(f"- `{provide}`")

                lines.append("")
                step_num += 1

        return "\n".join(lines)

    def describe_markdown(self, inputs: Dict[str, Any], viz_type: str = 'mermaid',
                        output_folder: Optional[Path] = None) -> str:
        """Generate markdown description including workflow visualization."""
        lines = [
            f"# Recipe: {self.name}",
            "",
            self.description,
            "",
            "## Inputs",
            ""
        ]

        for section, values in inputs.items():
            if values:  # Only show sections with values
                lines.append(f"### {section.title()}")
                lines.append("```yaml")
                for key, value in values.items():
                    lines.append(f"{key}: {value}")
                lines.append("```")
                lines.append("")

        lines.append("## Workflow Structure")
        lines.append(self.generate_mermaid_graph())
        lines.append("")

        lines.append("## Detailed Plan")
        lines.append(self.plan())

        return "\n".join(lines)

    def generate_mermaid_graph(self) -> str:
        """Generate Mermaid flowcharts showing workflow structure with separate diagrams for each group."""
        
        def _label_for(step) -> str:
            """Generate label for a step, escaping special characters."""
            lbl = getattr(step, "label", None) or step.name
            return str(lbl).replace('\n', '\\n').replace('"', "'")

        all_diagrams = []
        
        # Generate a diagram for each group
        for group_idx, group in enumerate(self.workflow.groups):
            group_name = group.name.replace('"', "'")

            lines = [f'### Workflow Group: {group_name}']

            [lines.append(chart_start) for chart_start in ['```mermaid', 'flowchart TD']]

            # Add group title as a comment
            if group.iterate_over:
                group_name += f" (iterates over {group.iterate_over})"

            # Create a title node for the group

            node_id = 0
            node_map = {}
            
            # Create nodes for steps in this group
            for step in group.steps:
                label = _label_for(step)
                node_map[label] = f's{node_id}'
                lines.append(f'  s{node_id}["{label}"]')
                node_id += 1
        
            # Add sequential flow within the group
            steps = list(group.steps)
            for i in range(len(steps) - 1):
                current_label = _label_for(steps[i])
                next_label = _label_for(steps[i + 1])
                current_id = node_map[current_label]
                next_id = node_map[next_label]
                lines.append(f'  {current_id} --> {next_id}')

            # Add dependency connections within the group
            for i, step in enumerate(steps):
                for j, prev_step in enumerate(steps):
                    if j < i and set(step.needs) & set(prev_step.provides):
                        step_label = _label_for(step)
                        prev_label = _label_for(prev_step)
                        step_id = node_map[step_label]
                        prev_id = node_map[prev_label]
                        lines.append(f'  {prev_id} -.-> {step_id}')



            lines.append('```')
            lines.append('')  # Empty line between diagrams

            all_diagrams.extend(lines)
    
        # Generate an overview diagram showing group relationships
        if len(self.workflow.groups) > 1:
            all_diagrams.extend(['## Workflow Overview', ''])
            all_diagrams.extend(['```mermaid', 'flowchart LR'])

            # Create nodes for each group
            for i, group in enumerate(self.workflow.groups):
                group_name = group.name.replace('"', "'")
                if group.iterate_over:
                    group_name += f"\\n(iterates over {group.iterate_over})"
                all_diagrams.append(f'  G{i}["{group_name}"]')

            # Connect groups in sequence
            for i in range(len(self.workflow.groups) - 1):
                all_diagrams.append(f'  G{i} --> G{i+1}')

            # Add styling for iteration groups in overview
            iter_groups = []
            for i, group in enumerate(self.workflow.groups):
                if group.iterate_over:
                    iter_groups.append(f'G{i}')


            all_diagrams.append('```')
            all_diagrams.append('')

        return '\n'.join(all_diagrams)