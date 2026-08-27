from collections import defaultdict


try:
    from PyQt6 import QtCore, QtGui, QtWidgets
    from PyQt6.QtCore import Qt
    QT_VERSION = 6
except ImportError:
    try:
        from PyQt5 import QtCore, QtGui, QtWidgets
        from PyQt5.QtCore import Qt
        QT_VERSION = 5
    except ImportError:
        QT_VERSION = None
        QtCore = None
        QtGui = None
        QtWidgets = None
        Qt = None
        

## Cursor type mode
if QT_VERSION == 5:
    pointerType = defaultdict(lambda: Qt.ArrowCursor)
    pointerType['CUTTING_EMPTY'] = Qt.ForbiddenCursor
    pointerType['CUTTING_INITIAL'] = Qt.CrossCursor
    pointerType['CUTTING_WIDTH'] = Qt.CrossCursor
    pointerType['CUTTING_DIRECTION'] = Qt.CrossCursor
    pointerType['CUTTING_MOVE'] = Qt.OpenHandCursor
    pointerType['RESOLUTION'] = Qt.BlankCursor
elif QT_VERSION == 6:
    pointerType = defaultdict(lambda: Qt.CursorShape.ArrowCursor)
    pointerType['CUTTING_EMPTY'] = Qt.CursorShape.ForbiddenCursor
    pointerType['CUTTING_INITIAL'] = Qt.CursorShape.CrossCursor
    pointerType['CUTTING_WIDTH'] = Qt.CursorShape.CrossCursor
    pointerType['CUTTING_DIRECTION'] = Qt.CursorShape.CrossCursor
    pointerType['CUTTING_MOVE'] = Qt.CursorShape.OpenHandCursor
    pointerType['RESOLUTION'] = Qt.CursorShape.BlankCursor
else:
    pointerType = defaultdict(lambda: None)