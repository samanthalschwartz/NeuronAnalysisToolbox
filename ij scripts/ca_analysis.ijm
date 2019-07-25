run("FeatureJ Derivatives", "x-order=0 y-order=0 z-order=2 smoothing=2.0");
run("FeatureJ Laplacian", "compute smoothing=1.0");
run("Macro...", "code=v=-v stack");