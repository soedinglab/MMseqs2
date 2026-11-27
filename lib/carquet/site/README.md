# carquet website

Static GitHub Pages site (no build step): `index.html`, `style.css`, `app.js`, `assets/`. Deployed by `.github/workflows/pages.yml` on every push to `main` that touches `site/` (enable once in Settings → Pages → Source: GitHub Actions).

The benchmark charts and headline numbers are computed in the browser from `bench.json`, an unmodified `run_benchmark.py` report. To refresh them, rerun the benchmark and copy the new report over it:

```bash
python3 benchmark/run_benchmark.py --skip-xlarge --fastparquet --polars
cp benchmark/bench_<version>_<os>_<arch>_<date>.json site/bench.json
```

Local preview: `python3 -m http.server -d site 8765`, then open http://localhost:8765 (opening `index.html` from disk cannot fetch `bench.json`).
