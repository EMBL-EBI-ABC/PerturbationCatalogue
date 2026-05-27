from datetime import datetime
from pathlib import Path
from threading import Lock


LOG_SEPARATOR_WIDTH = 80
PRINT_LOCK = Lock()


def _ensure_log_file(log_file: str | Path) -> Path:
    """Normalize the log file path and ensure its parent directory exists."""
    log_file = Path(log_file).resolve()
    log_file.parent.mkdir(parents=True, exist_ok=True)
    return log_file


def append_log_block(log_file: str | Path, title: str, *lines: str) -> None:
    """Append a formatted multi-line status block to a log file."""
    timestamp = datetime.now().isoformat(timespec="seconds")
    separator = "=" * LOG_SEPARATOR_WIDTH
    log_file = _ensure_log_file(log_file)
    with log_file.open("a", encoding="utf-8") as handle:
        handle.write(f"\n[{timestamp}] {separator}\n")
        handle.write(f"[{timestamp}] {title}\n")
        for line in lines:
            handle.write(f"[{timestamp}] {line}\n")
        handle.write(f"[{timestamp}] {separator}\n")


def append_log_line(log_file: str | Path, message: str) -> None:
    """Append one timestamped message line to a log file."""
    timestamp = datetime.now().isoformat(timespec="seconds")
    log_file = _ensure_log_file(log_file)
    with log_file.open("a", encoding="utf-8") as handle:
        handle.write(f"[{timestamp}] {message}\n")


def print_status_block(
    log_file: str | Path,
    title: str,
    *lines: str,
    echo: bool = True,
) -> None:
    """Write a status block to a log file and optionally echo it to stdout."""
    separator = "=" * LOG_SEPARATOR_WIDTH
    with PRINT_LOCK:
        append_log_block(log_file, title, *lines)
        if not echo:
            return
        print(f"\n{separator}")
        print(title)
        for line in lines:
            print(line)
        print(separator)