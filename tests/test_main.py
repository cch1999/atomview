"""Tests for the atomview package entry point."""

from unittest.mock import patch

from atomview import main


def test_main() -> None:
    """Test the main function outputs the expected message."""
    with patch("builtins.print") as mock_print:
        main()
        mock_print.assert_called_once_with("Hello from atomview!")
