from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, List, Optional, Callable
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
        for part in parts[1:]:
            value = value[part]

        # Get actual parameter value
        if isinstance(value, dict) and param_name in value:
            value = value[param_name]

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
    executor: Optional[FunctionExecutor] = None

    @classmethod
    def from_yaml(cls, name: str, data: Dict[str, Any]) -> 'Step':
        """Create a Step instance from YAML data."""
        params = {
            param_name: InputField(**param_def)
            for param_name, param_def in data.get('params', {}).items()
        }

        executor = None
        if 'algorithm' in data and 'function' in data['algorithm']:
            executor = FunctionExecutor(
                function_path=data['algorithm']['function'],
                function_args=data['algorithm'].get('args', []),
                params=params
            )

        return cls(
            name=name,
            description=data.get('description', ''),
            params=params,
            needs=data.get('needs', []),
            provides=data.get('provides', []),
            executor=executor
        )

    def execute(self, state: Dict[str, Any]) -> None:
        """Execute the step and update state."""
        if not self.executor:
            raise ValueError(f"Step {self.name} has no executor")

        try:
            result = self.executor.execute(state)
        except Exception as e:
            logger.error(f"Failed to execute step {self.name}: {str(e)}")
            raise

        # Update state with results
        if isinstance(result, tuple) and len(result) == len(self.provides):
            for val, output_name in zip(result, self.provides):
                state[output_name] = val
        elif len(self.provides) == 1:
            state[self.provides[0]] = result
        else:
            raise ValueError(f"Unexpected result format from {self.name}")


@dataclass
class Recipe:
    """Main recipe class that coordinates steps and manages execution."""
    name: str
    description: str
    steps: List[Step]
    input_schema: Dict[str, Dict[str, InputField]]
    outputs: List[str]
    
    def __post_init__(self):
        """Create dict of steps for faster lookup."""
        self.steps_dict = {step.name: step for step in self.steps}

    @classmethod
    def load(cls, recipe_name: str, base_path: Path) -> 'Recipe':
        """
        Loads a recipe definition, its input schema, and operations based on the given
        recipe name and base path, returning an initialized instance of the Recipe class.

        This method reads and processes YAML files corresponding to the recipe definition,
        schema, and operations. It constructs a Recipe instance using this data.

        Args:
            recipe_name: The name of the recipe to load.
            base_path: The base directory path where recipe, schema, and operations YAML files
                are located.

        Returns:
            Recipe: An instance of the Recipe class constructed from the loaded recipe data.
        """
        # Load recipe definition
        with open(base_path / "recipes" / f"{recipe_name}.yml") as f:
            recipe_data = yaml.safe_load(f)

        # Load input schema
        with open(base_path / "schema" / f"{recipe_name}_schema.yml") as f:
            schema_data = yaml.safe_load(f)['schema']
            input_schema = {
                section: {
                    field_name: InputField(**field_def)
                    for field_name, field_def in fields.items()
                }
                for section, fields in schema_data.items()
            }

        # Load operations
        with open(base_path / "schema" /  f"{recipe_name}_operations.yml") as f:
            operations_data = yaml.safe_load(f)['schema']['operations']
            steps = [
                Step.from_yaml(step_name, operations_data[step_name])
                for step_name in recipe_data['steps']
            ]

        return cls(
            name=recipe_name,
            description=recipe_data.get('description', ''),
            steps=steps,
            input_schema=input_schema,
            outputs=recipe_data.get('outputs', [])
        )

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

    def get_param_value(self, param_def, state):
        """
        Retrieves a parameter value from the provided state based on the definition of the parameter.

        This function determines the value of a parameter from a given state object. The retrieval logic
        is based on where the parameter is defined to exist (e.g., within the `state` dictionary or its
        nested structures). It offers backward compatibility for older parameter definitions.

        Args:
            param_def (dict): The definition of the parameter, including information
                about its source location within the state.
            state (dict): The state object containing current values of parameters
                and possibly nested inputs used to resolve the parameter value.

        Returns:
            Any: The resolved value of the parameter based on its definition in the
                 state.
        """
        if param_def.get('from') == 'state':
            return state[param_def['name']]
        elif param_def.get('from', '').startswith('inputs.'):
            section, key = param_def['from'].split('.')[1:]
            return state['inputs'][section][key]
        else:
            # Default to old behavior for backward compatibility
            return state.get(param_def['name'])

    def validate_state(self, state: Dict[str, Any], step: Step) -> List[str]:
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

    def run(self, inputs: Dict[str, Any], outdir: Path) -> Dict[str, Any]:
        """
        Executes a sequence of steps in the recipe pipeline, performing state validation,
        step execution, and output verification at each step, and ultimately produces the
        final outputs of the recipe.

        Args:
            inputs (Dict[str, Any]): Input parameters or data required by the recipe.
            outdir (Path): Directory path where output data of the recipe will be stored.

        Returns:
            Dict[str, Any]: A dictionary containing the final output values of the recipe.

        Raises:
            ValueError: If state validation fails before executing a step, or if output
                verification fails after executing a step or upon finalizing the recipe.
            Exception: If an exception occurs during step execution, the error is logged
                and re-raised.
        """
        state = {
            'outdir': outdir,
            'inputs': inputs
        }

        logger.info(f"Starting recipe: {self.name}")

        for step in self.steps:
            logger.info(f"Executing step: {step.name}")
            
            # Validate state before execution
            if errors := self.validate_state(state, step):
                raise ValueError(f"Step {step.name} state validation failed:\n" + "\n".join(errors))

            # Execute step
            try:
                step.execute(state)
            except Exception as e:
                logger.error(f"Failed to execute step {self.name}: {str(e)}")
                raise

            # Verify step outputs
            if errors := self.verify_outputs(state, step):
                raise ValueError(f"Step {step.name} output verification failed:\n" + "\n".join(errors))

            logger.info(f"Completed step: {step.name}")

        # Verify recipe outputs and prepare return value
        if errors := self.verify_recipe_outputs(state):
            raise ValueError("\n".join(errors))

        recipe_outputs = {out: state[out] for out in self.outputs}
        logger.info(f"Completed recipe: {self.name}")

        return recipe_outputs

    def plan(self) -> str:
        """
        Generates a formatted plan as a string representation from the provided steps. It iterates through
        each defined step, appending its name, description, parameters, and any provided outputs to a
        formatted string. Each step is properly numbered and separated for better readability.

        Returns:
            str: A formatted string combining all steps, detailing their relevant information in a structured format.
        """
        lines = []

        for i, step in enumerate(self.steps, 1):
            lines.append(f"### {i}. {step.name}")

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

            lines.append("")  # Empty line between steps

        return "\n".join(lines)

    def visualize_dependencies_networkx(self, output_folder: Optional[Path] = None) -> str:
        """
        Visualizes the dependencies of steps in the current process using NetworkX and Matplotlib.
        Generates a directed graph where nodes represent steps, and edges indicate dependencies
        between steps based on provided and required resources.

        Args:
            output_folder (Optional[Path]): A folder where the generated dependency graph image
                will be saved. If not provided, the method will not save the visualization.

        Returns:
            str: A markdown string referencing the saved image if an output folder is provided,
                or a message indicating that no output folder was given.
        """
        import networkx as nx
        import matplotlib.pyplot as plt

        G = nx.DiGraph()

        # Add nodes
        for step in self.steps:
            G.add_node(step.name)

        # Add edges
        for step in self.steps:
            for prev_step in self.steps:
                if set(step.needs) & set(prev_step.provides):
                    G.add_edge(prev_step.name, step.name)

        # Create the plot
        plt.figure(figsize=(12, 8))
        pos = nx.spring_layout(G, k=1, iterations=50)
        nx.draw(G, pos, with_labels=True, node_color='lightblue',
                node_size=2000, font_size=10, font_weight='bold',
                arrows=True, edge_color='gray')
        plt.title(f"Dependencies for {self.name}")

        # Save if output folder is provided
        if output_folder:
            output_folder = Path(output_folder)
            figure_path = output_folder / f"{self.name}_dependencies.png"
            plt.savefig(figure_path, bbox_inches='tight', dpi=300)
            plt.close()

            # Return markdown that references the saved image
            return f"![{self.name} Dependencies]({figure_path.name})"
        else:
            plt.close()
            # Return a message if no output folder was provided
            return "*NetworkX visualization requires an output folder to save the figure.*"

    def generate_mermaid_graph(self) -> str:
        """
        Generates a Mermaid graph representation based on the defined steps and their
        dependencies.

        This function constructs a flowchart in the Mermaid graph syntax by creating
        nodes and defining connections based on the main flow and dependency
        relationships between the steps.

        Returns:
            str: The Mermaid graph representation as a string.
        """
        lines = ['```mermaid', 'graph TD;']

        # Add nodes
        for step in self.steps:
            lines.append(f'    {step.name}[{step.name}];')

        # Add main flow connections
        for i in range(len(self.steps) - 1):
            lines.append(f'    {self.steps[i].name} --> {self.steps[i+1].name};')

        # Add dependency connections
        for i, step in enumerate(self.steps):
            for j, prev_step in enumerate(self.steps):
                if j < i and set(step.needs) & set(prev_step.provides):
                    lines.append(f'    {prev_step.name} -.-> {step.name};')

        lines.append('```')
        return '\n'.join(lines)

    def describe_markdown(self, inputs: Dict[str, Any], viz_type: str = 'mermaid',
                          output_folder: Optional[Path] = None) -> str:
        """
        Generates a markdown description of the recipe, including its inputs, execution plan, and
        dependency graph visualized using the specified method.

        Args:
            inputs (Dict[str, Any]): A dictionary of inputs categorized into sections where keys are
                section names, and values are further dictionaries mapping input names to their values.
            viz_type (str): The type of visualization for the dependency graph. Supported options are
                'mermaid', 'ascii', and 'networkx'. Default is 'mermaid'.
            output_folder (Optional[Path]): The folder path where graph visualizations will be saved
                (if required by the selected visualization method).

        Returns:
            str: The generated markdown description as a string.

        Raises:
            ValueError: If an unsupported visualization type is provided in `viz_type`.
        """
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
                lines.append("```yaml")  # Using yaml for better formatting
                for key, value in values.items():
                    lines.append(f"{key}: {value}")
                lines.append("```")
                lines.append("")

        lines.append("## Execution Plan")
        lines.append(self.plan())

        lines.append("## Dependency Graph")

        # Choose visualization type
        viz_methods = {
            'mermaid': lambda: self.generate_mermaid_graph(),
            'ascii': lambda: self.generate_ascii_graph(),
            'networkx': lambda: self.visualize_dependencies_networkx(output_folder)
        }

        if viz_type not in viz_methods:
            raise ValueError(f"Unsupported visualization type. Choose from: {', '.join(viz_methods.keys())}")

        lines.append(viz_methods[viz_type]())

        return "\n".join(lines)

    def generate_ascii_graph(self) -> str:
        """
        Generates an ASCII representation of the dependency graph for steps.

        This function creates a visual representation of how the steps in a particular
        workflow are dependent on one another. Each step is visualized as part of a tree,
        showing parent-child relationships between dependency nodes. The root nodes,
        representing steps with no dependencies, are identified and used to construct
        the graph recursively.

        Returns:
            str: A string containing the ASCII representation of the dependency graph.
        """
        lines = []
        indent = "  "

        def find_dependencies(step_name: str, depth: int = 0) -> None:
            lines.append(f"{indent * depth}└─ {step_name}")
            for next_step in self.steps:
                if any(dep in self.steps_dict[step_name].provides
                      for dep in next_step.needs):
                    find_dependencies(next_step.name, depth + 1)

        # Find root nodes (no dependencies)
        roots = [step.name for step in self.steps if not step.needs]
        for root in roots:
            find_dependencies(root)

        return "\n".join(lines)