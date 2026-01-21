"""Toy tests for package-template."""

from unittest.mock import patch

import pytest

from package_template import main
from package_template.cli import main as cli_main


def test_main():
    """Test the main function outputs the expected message."""
    with patch("builtins.print") as mock_print:
        main()
        mock_print.assert_called_once_with("Hello from package-template!")


def test_cli_default(capsys):
    """Test CLI with default arguments."""
    with patch("sys.argv", ["cli"]):
        cli_main()
    captured = capsys.readouterr()
    assert "Hello, World!" in captured.out


def test_cli_custom_name(capsys):
    """Test CLI with custom name argument."""
    with patch("sys.argv", ["cli", "--name", "Alice"]):
        cli_main()
    captured = capsys.readouterr()
    assert "Hello, Alice!" in captured.out


def test_cli_help(capsys):
    """Test CLI help message."""
    with pytest.raises(SystemExit):
        with patch("sys.argv", ["cli", "--help"]):
            cli_main()
    # Help should be printed to stdout
    captured = capsys.readouterr()
    assert "Example CLI program" in captured.out or "package-template" in captured.out
