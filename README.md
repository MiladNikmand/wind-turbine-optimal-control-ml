# wind-turbine-optimal-control-ml
Master's thesis project on wind turbine MPPT control. Wind speed is first predicted using machine learning, then a proposed DDP-HJB solver computes the optimal control offline over that predicted horizon — so the optimal solution is ready and applied in real time, enabling near-optimal performance without heavy onboard computation.
