[//]: # (    Note that this file is written in Markdown, and is better viewed with a Markdown viewer/editor, such as https://dillinger.io     )



# Time-Dependent Channel Flows - D2Q9 Lattice Boltzmann implementation

This repository contains code to run Poiseuille-flow lattice Boltzmann simulations in a browser, with a visual interface.

Implementation written in JavaScript (ES6), with [Three.js](https://threejs.org/) for the real-time render and vector maths, and [Chart.js](https://www.chartjs.org/) for the runtime-generated charts.

## Getting started

Cloning the repository locally and opening [index.html](index.html) in a browser (any modern browser should work; tested on Firefox 93.0 and Chrome 94.0 on Windows 10 (64-bit)), either using a `file:///` protocol or with a local server on `localhost`, will bring up a 100x100 pulsatile Poiseuille flow D2Q9 simulation. The simulation parameters can be tweaked in `index.html` directly.

To set up a separate simulation, the following files must be included in the HTML document `<head>`: Three.js version r128 (from CDN here), Chart.js 3.5.1 (from CDN), `plots.js` (only if plotting capabilities are needed), `renderer.js` (only if needing to render the simulation in real-time), `ui.js` (utility script to add buttons), and the main D2Q9 simulation implementation via `lbm.js`.

```html
<head>
    ...
    <script src='https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js'></script>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@3.5.1/dist/chart.min.js"></script>
    <script src='plots.js'></script>
    <script src='renderer.js'></script>
    <script src='ui.js'></script>
    <script src='lbm.js'></script>
    ...
</head>
```

To initialize a simulation, an instance of `D2Q9.Parameters` can be created and filled:

```js
let simParams = new D2Q9.Parameters();
simParams.width = 20; // number of lattice sites horizontally
simParams.height = 20; // number of lattice sites across the channel
simParams.kinematicViscosity = 0.5; // kinematic velocity nu
simParams.initialDensity = 1.0; // initial density for all sites rho
simParams.force = (time) => new THREE.Vector2(5e-8, 0); // optional time-dependent force
```

Then, it can be passed to the constructor of `D2Q9`:

```js
let sim = new D2Q9(simParams);
```

Boundaries (with half-way bounce-back boundary condition) can be setup on the `D2Q9` instance with `addHorizontalBoundary` and `addCircularBoundary`:

```js
sim.addHorizontalBoundary(0); // bottom boundary
sim.addHorizontalBoundary(simParams.height - 1); // top boundary
```

Any edges without any boundaries will be periodically connected to the other side of the lattice on both axes.

Finally, the simulation can be progressed by a single time-step with `step()`:

```js
// Progress the simulation by 100 time-steps
for(let i = 0; i < 100; ++i) {
    sim.step();
}
```
