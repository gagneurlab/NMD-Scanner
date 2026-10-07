"""polars-bio, imported without the change it makes to the root logger."""

import logging

# Importing polars-bio calls logging.basicConfig() and sets the root level to WARNING. Both are
# undone here, so that importing nmd_scanner leaves the logging of the caller as it was.
_root = logging.getLogger()
_handlers = _root.handlers[:]
_level = _root.level
import polars_bio as pb  # noqa: E402

_root.handlers[:] = _handlers
_root.setLevel(_level)

__all__ = ["pb"]
