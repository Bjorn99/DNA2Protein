# Contributing Guidelines for DNA2PROTEIN

## Welcome Contributors

Thank you for your interest in contributing to DNA2PROTEIN. This document provides guidelines for contributing to this project.

## Table of Contents
1. [Code of Conduct](#code-of-conduct)
2. [Getting Started](#getting-started)
3. [How to Contribute](#how-to-contribute)
4. [Development Setup](#development-setup)
5. [Pull Request Process](#pull-request-process)
6. [Style Guidelines](#style-guidelines)
7. [Testing Guidelines](#testing-guidelines)
8. [Documentation](#documentation)

## Code of Conduct

By participating in this project, you agree to maintain a respectful and collaborative environment. We expect all contributors to:
- Use welcoming language
- Be respectful of differing viewpoints
- Accept constructive criticism gracefully
- Focus on what is best for the community
- Show empathy towards other community members

## Getting Started

1. Fork the repository
2. Clone your fork:
```bash
git clone https://github.com/your-username/DNA2PROTEIN.git
cd DNA2PROTEIN
```
3. Set up your development environment following the [Development Setup](#development-setup) section

## How to Contribute

### Reporting Bugs
- Use the GitHub Issues section
- Check if the bug has already been reported
- Include detailed steps to reproduce the bug
- Provide system information and relevant details

### Suggesting Enhancements
- First, check if the enhancement has already been suggested
- Use the GitHub Issues section with the "enhancement" label
- Provide clear and detailed explanation of the proposed feature
- Include any relevant scientific references or documentation

### Code Contributions

Areas where contributions are particularly welcome:
- Improving signal peptide prediction algorithms
- Adding new sequence analysis features
- Enhancing visualization capabilities
- Optimizing performance
- Improving documentation
- Adding test cases

## Development Setup

1. Ensure you have Python 3.11+ installed
2. Install poetry for dependency management:
```bash
pip install poetry
```
3. Install dependencies:
```bash
poetry install
```
4. Set up pre-commit hooks:
```bash
pre-commit install
```

## Pull Request Process

1. Create a new branch for your feature:
```bash
git checkout -b feature/your-feature-name
```

2. Make your changes, following our [Style Guidelines](#style-guidelines)

3. Write or update tests as needed

4. Update documentation as needed

5. Commit your changes:
```bash
git add .
git commit -m "Description of changes"
```

6. Push to your fork:
```bash
git push origin feature/your-feature-name
```

7. Create a Pull Request through GitHub

## Style Guidelines

### Python Code Style
- Follow PEP 8 guidelines
- Use type hints where applicable
- Maximum line length: 88 characters
- Use descriptive variable names
- Document functions using docstrings

### JavaScript Style
- Use ES6+ features where appropriate
- Follow Airbnb JavaScript Style Guide
- Use meaningful variable and function names

### HTML/CSS Style
- Follow BEM methodology for CSS
- Use semantic HTML elements
- Maintain consistent indentation

## Testing Guidelines

1. Write tests for new features
2. Ensure all tests pass before submitting PR:
```bash
poetry run pytest
```
3. Include both unit tests and integration tests where applicable
4. Aim for meaningful test coverage

## Documentation

1. Update docstrings for any new functions
2. Update README.md if adding new features
3. Include inline comments for complex logic
4. Update requirements if adding new dependencies

## Additional Notes

### Commit Messages
- Use clear and descriptive commit messages
- Reference issues and pull requests liberally
- Format: `[type]: Description of changes`
  - Types: feat, fix, docs, style, refactor, test, chore

### Version Control
- Keep pull requests focused on a single feature/fix
- Rebase your branch before submitting PR
- Squash commits if necessary

### Questions?
Feel free to create an issue or contact the maintainers if you have any questions.

## License
By contributing to DNA2PROTEIN, you agree that your contributions will be licensed under the project's license.

Thank you for helping improve DNA2PROTEIN!
