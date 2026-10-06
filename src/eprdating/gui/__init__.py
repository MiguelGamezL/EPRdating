"""Interactive interface for Jupyter and Google Colab.

    from eprdating.gui import app
    app()

Three tabs: the equivalent dose from EPR spectra, the sediment U, Th and K
from HPGe spectra, and the age. A fitted De and an analysed sediment flow
into the Age tab. Needs ``ipywidgets`` and ``matplotlib``
(``pip install "eprdating[gui]"``; both come with Google Colab).

Each tab is also available alone (:class:`DePanel`, :class:`GammaPanel`,
:class:`AgePanel`) and can be driven from code, which is how it is tested.
"""

from __future__ import annotations

try:
    import ipywidgets as _w
except ImportError as e:  # pragma: no cover
    raise ImportError('the interface needs ipywidgets: pip install "eprdating[gui]"') from e

from .age import AgePanel
from .de import DePanel
from .gamma import GammaPanel


class App:
    """The three panels in tabs. ``app.de``, ``app.gamma`` and ``app.age``
    are the panels; ``app.widget`` is what is displayed."""

    def __init__(self) -> None:
        self.age = AgePanel()
        self.de = DePanel(on_de=self.age.set_De)
        self.gamma = GammaPanel(on_sediment=self.age.set_sediment)
        self.widget = _w.Tab(children=[self.de.widget, self.gamma.widget, self.age.widget])
        for i, title in enumerate(("1 · Equivalent dose", "2 · Sediment", "3 · Age")):
            self.widget.set_title(i, title)
        from .. import __version__

        self.header = _w.HTML(f"<b>EPRdating {__version__}</b> &nbsp; "
                              "<small>ESR dating of tooth enamel · "
                              "<a href='https://github.com/MiguelGamezL/EPRdating' target='_blank'>documentation</a>"
                              "</small>")
        self.box = _w.VBox([self.header, self.widget])

    def _ipython_display_(self):
        from IPython.display import display

        display(self.box)


def app() -> App:
    """Open the interface (in a notebook cell: ``app()``)."""
    return App()


__all__ = ["AgePanel", "App", "DePanel", "GammaPanel", "app"]
