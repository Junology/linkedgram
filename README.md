# linkedgram
A program written in Haskell to draw link projection diagrams.

## License
The project is released under BSD 3 License; see LICENSE.

## Build
Since the project uses Gtk2Hs (https://github.com/gtk2hs/gtk2hs) for GUI, you may need to install it and related libraries before building this repository.
For this, follow the instruction on their project page or HaskellWiki (https://wiki.haskell.org/Gtk2Hs/Installation).
MacOS users may need to see https://wiki.haskell.org/Gtk2Hs/Mac.

Now, the easiest way to build the program is to use Stack (https://www.haskellstack.org).
So just type the folloing on your terminal:
```
stack install
```

## Features

- Draw link projection diagrams in Draw mode (see Warning below).
- When you finish drawing (by Right-click or toggling "Draw mode" button), crossings are automatically detected.
- Double-click on crossings for crossing-change.
- Right-click on crossing for smoothing.
- Exporting to LaTeX is supported using TikZ package.

### Warning
If the following materials are overlapped, then the diagram is "ill-formed;" i.e. it is no longer valid.

- Control points you put in Draw mode.
- Convex hull of control points involved with a crossing.

To avoid this bad situation, you need to put enough numbers of control points.
