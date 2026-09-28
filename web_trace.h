#pragma once

#include <vector>

#include "simplify_geometry.h"

// Web-trace emission: mirrors the stab loop of simplify_core.h but records
// every intermediate value (P, Gi, F, S) for the web visualizer. With
// --json-stream it emits NDJSON (header, one prefix per line, done);
// otherwise it emits a single JSON object. Defined in web_trace.cpp.
std::vector<Point> simplify_web(const std::vector<Point>& stream,
                                double EPSILON, double DELTA);
