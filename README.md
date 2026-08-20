# Roman_GRS_footprint

This repository contains tools to simulate, use, and analyze the Roman GRS footprint.

## Documentation

The Sphinx documentation is in `docs/` and is configured to document the `rstgrs_footprint` package.

## Development

Install the package in editable mode with:

```bash
pip install -e .
```

To install documentation dependencies:

```bash
pip install -e ".[docs]"
```

## Optional Dependencies

### roman_gdps_optical_model

Some functionality in this package requires the `roman_gdps_optical_model` module, which is available in a private GitHub repository. To use features that depend on this module:

1. Ensure you have access to the private repository
2. Install it using:
   ```bash
   pip install git+https://github.com/roman-grs-pit/roman_gdps_optical_model.git
   ```
   (Requires appropriate GitHub credentials or authentication token)

3. Set the required environment variable:
   ```bash
   export ROMAN_GDPS_OPTICAL_MODEL_CONFIG=/path/to/config/file
   ```

Note: The standard test suite does not require this dependency. It is only needed for specific functionality related to optical modeling and footprint analysis.
