try:
    from PyQt6 import QtCore, QtGui, QtWidgets
    from PyQt6.QtCore import Qt
    QT_VERSION = 6
except ImportError:
    try:
        from PyQt5 import QtCore, QtGui, QtWidgets
        from PyQt5.QtCore import Qt
        QT_VERSION = 5
    except ImportError as exc:
        raise ImportError(
            "MJOLNIR requires either PyQt5 or PyQt6. "
            "Install one with `pip install MJOLNIR[qt5]` "
            "or `pip install MJOLNIR[qt6]`."
        ) from exc