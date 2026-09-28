# Interactive temporal–comorbidity architecture associated with low grip strength

This repository contains the interactive three-dimensional visualization of the integrated temporal–comorbidity architecture associated with low grip strength.

The visualization accompanies the analyses presented in Fig. 4 of the manuscript and is intended as an interactive extension of the static figure rather than an independent analysis.

## Interactive visualization

The browser-based visualization is available at:

https://leonalinmeixiu-dot.github.io/LowGrip-Temporal-Network/

Users can rotate, zoom and explore individual diseases and their relationships within the integrated temporal–comorbidity architecture.

## Repository contents

- `index.html` — browser-based interactive 3D visualization
- `nodes.json` — processed node-level information used by the visualization
- `edges.json` — processed edge-level information used by the visualization
- `code/` — scripts used to generate the processed network files and visualization inputs

## Visualization framework

Nodes represent diseases included in the multimorbidity and temporal analyses.

Node colors indicate disease categories or network communities as defined in the study.

The spatial organization integrates disease-community structure with temporal ordering relative to low grip strength, including upstream, downstream and associated diseases.

The visualization does not introduce additional statistical analyses beyond those reported in the manuscript.

## Data availability note

Only processed network-level information required for visualization is included in this repository. Individual-level participant data are not included.

## Requirements

No software installation is required to explore the interactive visualization. A modern web browser with WebGL support is recommended.

## Citation

If you use this resource, please cite the accompanying manuscript. Full citation details will be added upon publication.

## License

The code and visualization components in this repository are released under the MIT License. See the `LICENSE` file for details.
