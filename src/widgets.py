import asyncio
import dataclasses
from typing import Any, Sequence, Iterable

import numbers
import ipywidgets as widgets
from ipysheet import Cell, Sheet


def determine_cell_type(value) -> str:
    if isinstance(value, bool):
        return 'checkbox'
    elif isinstance(value, numbers.Number):
        return 'numeric'
    elif isinstance(value, widgets.Widget):
        return 'widget'
    return 'text'


def base_field_transformer(dc: dataclasses.dataclass, field_name: str) -> dict[str, Any]:
    return dict(value=getattr(dc, field_name), type=determine_cell_type(dc))


def sheet_from_dataclass(
        data: Sequence[dataclasses.dataclass],
        *,
        fields: list[str] | None = None,
        field_transformers: dict[str, callable] | None = None,
        cell_kwargs: dict | None = None,
        sheet_kwargs: dict | None = None,
) -> Sheet:
    if not data:
        return Sheet()

    fields = fields if fields is not None else [
        field.name for field in dataclasses.fields(data[0])
    ]

    cell_kwargs = cell_kwargs or {}
    cell_kwargs.setdefault('read_only', True)
    sheet_kwargs = sheet_kwargs or {}
    field_transformers = field_transformers or {}

    def row_cells(dc: dataclasses.dataclass, row_idx: int) -> Iterable[Cell]:
        fields_by_name = {field.name: field for field in dataclasses.fields(dc)}
        dc_fields = [fields_by_name[name] for name in fields]
        return [
            Cell(
                row_start=row_idx,
                row_end=row_idx,
                column_start=col_idx,
                column_end=col_idx,
                squeeze_row=True,
                squeeze_column=True,
                **{
                    **cell_kwargs,
                    **base_field_transformer(dc, field.name),
                    **field_transformers.get(field.name, lambda *_args: {})(dc, field.name),
                }
            )
            for col_idx, field in enumerate(dc_fields)
        ]

    cells = [cell for row_id, dc in enumerate(data) for cell in row_cells(dc, row_id)]

    return Sheet(
        rows=len(data),
        row_headers=True,
        columns=len(fields),
        column_headers=fields,
        cells=cells,
        **sheet_kwargs,
    )


#
#  ipywidgets debouncing
#  (https://ipywidgets.readthedocs.io/en/8.1.5/examples/Widget%20Events.html#debouncing)
#

class TimerAsync:
    def __init__(self, timeout, callback):
        self._timeout = timeout
        self._callback = callback
        self._task = None

    async def _job(self):
        await asyncio.sleep(self._timeout)
        self._callback()

    def start(self):
        self._task = asyncio.ensure_future(self._job())

    def cancel(self):
        self._task.cancel()


def debounce(wait, *, timer_factory=TimerAsync):
    """ Decorator that will postpone a function's
        execution until after `wait` seconds
        have elapsed since the last time it was invoked. """

    def decorator(fn):
        timer = None

        def debounced(*args, **kwargs):
            nonlocal timer

            def call_it():
                fn(*args, **kwargs)

            if timer is not None:
                timer.cancel()
            timer = timer_factory(wait, call_it)
            timer.start()

        return debounced

    return decorator
