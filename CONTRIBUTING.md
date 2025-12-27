# Contributing to LIPAC

We welcome contributions to LIPAC! This document provides guidelines for contributing to the project.

## How to Contribute

### Reporting Issues

1. Check if the issue has already been reported
2. Create a new issue with a clear title and description
3. Include:
   - LIPAC version
   - Python version
   - Operating system
   - Minimal reproducible example
   - Error messages (if any)

### Submitting Pull Requests

1. Fork the repository
2. Create a new branch for your feature (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests if applicable
5. Ensure all tests pass
6. Commit your changes (`git commit -m 'Add amazing feature'`)
7. Push to your branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

### Code Style

- Follow PEP 8 guidelines
- Use meaningful variable and function names
- Add docstrings to all functions and classes
- Keep line length under 100 characters
- Use type hints where appropriate

### Testing

- Write tests for new functionality
- Ensure existing tests pass
- Aim for high test coverage
- Test with multiple Python versions (3.10+)

### Documentation

- Update documentation for new features
- Include docstrings in NumPy style
- Update README.md if needed
- Add examples for new functionality

## Development Setup

1. Clone the repository:
```bash
git clone https://github.com/takeshi-sato-dev/lipac.git
cd lipac
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install development dependencies:
```bash
pip install -e .[dev]
```

4. Run tests:
```bash
pytest
```

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers and help them get started
- Focus on constructive criticism
- Respect differing viewpoints and experiences

## Questions?

Feel free to open an issue for any questions about contributing.

Thank you for contributing to LIPAC!
