# Developer Guide: Creating New Recipes in the Recipe System

## Overview

The recipe system provides a flexible, YAML-based framework for defining data processing pipelines. A recipe consists of three main components:

1. `recipe.yml` - The main recipe definition
2. `recipe_operations.yml` - Operation/step definitions
3. `recipe_schema.yml` - Input validation rules

## File Structure and Relationships

### Base Recipe File (`recipe.yml`)

This is the main entry point that defines your recipe's structure


### Schema Definition (`recipe_schema.yml`)

Defines validation rules for all inputs


### Operations Definition (`recipe_operations.yml`)

Defines the implementation details for each step


## Creating a New Recipe

### Step 1: Plan Your Recipe

1. Identify the required inputs
2. Break down the process into discrete steps
3. Determine the dependencies between steps
4. Define the expected outputs

### Step 2: Create Recipe Files

1. Create `recipe.yml` with basic structure
2. Add input validation to `recipe_schema.yml`
3. Define operations in `recipe_operations.yml`

### Step 3: Implement Functions

Implement the Python functions referenced in your operations.

## Key Concepts

### Input Sources

- `inputs.cli`: Command-line arguments
- `inputs.dataset`: Dataset-specific inputs
- `inputs.config`: Configuration parameters
- `state`: Values from previous step outputs

### Supported Parameter Types

- Basic Types:
  - `str`: String values
  - `int`: Integer values
  - `float`: Floating-point values
  - `bool`: Boolean values
  
- Complex Types:
  - `path`: File system paths
  - `ndarray`: NumPy arrays
  - `dict`: Dictionaries
  
- List Types:
  - `list[type]`: Lists of specific types
  - Example: `list[str]`, `list[float]`

### Validation Rules

- Required Fields:
  - `required`: true/false
  - `description`: Documentation string
  - `type`: Data type

- Optional Fields:
  - `units`: Unit specification
  - `default`: Default value
  - `validate`: Validation rules

- Validation Options:
  - `min`/`max`: Numeric bounds
  - `min_length`/`max_length`: Length constraints
  - `pattern`: Regex patterns
  - `choices`: Allowed values
  - `required_keys`: Required dictionary keys

## Best Practices

### Recipe Organization

1. Group Related Steps
   - Order steps logically
   - Keep related operations together
   - Minimize cross-step dependencies

2. Input Management
   - Use clear, descriptive names
   - Group related inputs
   - Provide sensible defaults

3. Validation
   - Add comprehensive validation rules
   - Include units where applicable
   - Document constraints

### Step Design

1. Modularity
   - Keep steps focused and single-purpose
   - Make steps reusable when possible
   - Clear input/output interfaces

2. Dependencies
   - Minimize step dependencies
   - Document dependencies clearly
   - Use state management appropriately

3. Error Handling
   - Validate inputs thoroughly
   - Provide clear error messages
   - Handle edge cases

### Documentation

1. Descriptions
   - Clear step descriptions
   - Document parameter purposes
   - Explain validation rules

2. Examples
   - Provide usage examples
   - Document expected inputs/outputs
   - Include common patterns

## Testing and Debugging

1. Validate recipe structure
2. Test each step independently
3. Verify input validation rules
4. Check state management
5. Test full recipe execution

## Conclusion

The recipe system provides a flexible framework for creating data processing pipelines. 
By following these guidelines and best practices, you can create maintainable, 
well-documented recipes that are easy to understand and modify.
