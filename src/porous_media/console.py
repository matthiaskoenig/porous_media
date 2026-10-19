"""Rich console for logging."""

from rich.console import Console
from rich.theme import Theme

custom_theme = Theme(
    {
        "success": "green",
        "info": "blue",
        "warning": "orange3",
        "error": "red",
    }
)

console = Console(theme=custom_theme)
