"""
Centralized logging configuration for the orphan drug overlap pipeline.
"""

import logging
from pathlib import Path


def setup_logging(level: int = logging.INFO, 
                 format_str: str = "%(asctime)s %(levelname)s: %(message)s",
                 log_file: Path = None) -> None:
    """
    Configure logging for the entire pipeline.
    
    Args:
        level: Logging level (default: INFO)
        format_str: Log message format string
        log_file: Optional file to write logs to
    """
    handlers = [logging.StreamHandler()]
    
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=level,
        format=format_str,
        handlers=handlers,
        force=True  # Override any existing configuration
    )


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance for a specific module.
    
    Args:
        name: Logger name (typically __name__)
        
    Returns:
        Configured logger instance
    """
    return logging.getLogger(name) 