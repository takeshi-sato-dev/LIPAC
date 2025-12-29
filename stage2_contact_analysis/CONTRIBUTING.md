# Contributing to MIRAGE4

We welcome contributions to MIRAGE4! This document provides guidelines for contributing to the project.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/yourusername/mirage4.git
   cd mirage4
   ```
3. Create a development environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -r requirements.txt
   pip install -e .  # Install in development mode
   ```

## Development Process

1. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature-name
   ```

2. Make your changes and ensure:
   - Code follows PEP 8 style guidelines
   - All tests pass
   - Documentation is updated if needed

3. Run tests:
   ```bash
   pytest tests/
   ```

4. Check code style:
   ```bash
   flake8 mirage4/
   black --check mirage4/
   ```

5. Format code if needed:
   ```bash
   black mirage4/
   ```

## Submitting Changes

1. Commit your changes:
   ```bash
   git add .
   git commit -m "Descriptive commit message"
   ```

2. Push to your fork:
   ```bash
   git push origin feature-name
   ```

3. Submit a pull request through GitHub

## Pull Request Guidelines

- Provide a clear description of the changes
- Reference any related issues
- Include tests for new functionality
- Update documentation as needed
- Ensure all CI checks pass

## Code Style

- Follow PEP 8
- Use descriptive variable names
- Add docstrings to all functions and classes
- Keep functions focused and modular
- Maximum line length: 100 characters

## Testing

- Write tests for new features
- Maintain or improve code coverage
- Test edge cases and error conditions
- Use pytest for testing

## Documentation

- Update docstrings for API changes
- Update README.md if adding new features
- Add examples for new functionality
- Keep documentation clear and concise

## Reporting Issues

When reporting issues, please include:
- Python version
- MIRAGE4 version
- Minimal reproducible example
- Full error traceback
- Operating system

## Feature Requests

We welcome feature requests! Please:
- Check existing issues first
- Clearly describe the proposed feature
- Explain the use case
- Provide examples if possible

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Focus on constructive criticism
- Respect differing opinions

## Questions?

Feel free to open an issue for any questions about contributing!

## License

By contributing, you agree that your contributions will be licensed under the MIT License.