#!/usr/bin/env python3
"""
run_analysis.py — Orchestrator for the charm bulk workflow.

Reads a flow configuration YAML and runs the requested nodes:
  1. make_cutsets_cfgs  — Generate cutset YAML files
  2. mass_fit           — Fit invariant mass distributions (flarefly)
  3. cut_variation      — Corrected yields via minimisation
  4. data_driven_fraction — Prompt/FD fractions

Each node can take custom input paths; if not provided, auto-detects.

Usage:
    python3 run_analysis.py config.yml [--workers N] [--correlated]
"""
import os
import sys
import yaml
import argparse
import time
import concurrent.futures
import subprocess

# The python executable to use (with ROOT + flarefly)
PYTHON = os.environ.get("WF_PYTHON", "conda run -n alice python3")

script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(script_dir, "src")
sys.path.insert(0, src_dir)
from utils import logger, check_dir


def get_paths():
    """Map node names to their scripts."""
    return {
        "make_yaml": os.path.join(src_dir, "make_cutsets_cfgs.py"),
        "mass_fit": os.path.join(src_dir, "mass_fit.py"),
        "cut_variation": os.path.join(src_dir, "cut_variation.py"),
        "data_driven_fraction": os.path.join(src_dir, "data_driven_fraction.py"),
    }


def run_cmd(cmd, timeout=600):
    """Run a shell command with logging."""
    logger(f"Running: {cmd}", "COMMAND")
    result = subprocess.run(cmd, shell=True, timeout=timeout,
                            capture_output=True, text=True)
    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            logger(line, "OUTPUT")
    if result.returncode != 0:
        logger(f"Command failed (rc={result.returncode})", "ERROR")
        if result.stderr:
            for line in result.stderr.strip().split("\n")[-20:]:
                logger(f"  {line}", "ERROR")
        raise RuntimeError(f"Command failed: {cmd}")
    return result


def main():
    parser = argparse.ArgumentParser(description="Charm bulk workflow")
    parser.add_argument("flow_config", help="Flow configuration YAML file")
    parser.add_argument("--workers", "-w", type=int, default=1, help="Number of workers")
    parser.add_argument("--correlated", "-c", action="store_true", help="Correlated cut variation")
    parser.add_argument("--test", "-t", action="store_true", help="Test mode (3 pT bins)")
    args = parser.parse_args()

    start_time = time.time()

    # Load config
    with open(args.flow_config, "r") as f:
        config = yaml.safe_load(f)

    operations = config.get("operations", {})
    base_outdir = config.get("outdir", ".")
    suffix = config.get("suffix", "default")

    # Determine output directory
    corr_str = "_correlated" if args.correlated else "_combined"
    if "_correlated" not in base_outdir and "_combined" not in base_outdir:
        outdir = os.path.join(base_outdir, f"cutvar_{suffix}{corr_str}")
    else:
        outdir = base_outdir

    os.makedirs(outdir, exist_ok=True)

    # Copy config
    cfg_copy_dir = os.path.join(outdir, "config_flow")
    os.makedirs(cfg_copy_dir, exist_ok=True)
    cfg_copy = os.path.join(cfg_copy_dir, os.path.basename(args.flow_config))
    os.system(f"cp {args.flow_config} {cfg_copy}")

    # For testing: modify to use only first 3 pT bins
    if args.test:
        ptbins = config.get("ptbins", [])
        if len(ptbins) > 4:
            logger(f"TEST MODE: reducing from {len(ptbins)-1} to 3 pT bins", "WARNING")
            config["ptbins"] = ptbins[:4]  # Keep first 3 bins (4 edges)
            # Write modified config
            test_cfg_path = os.path.join(outdir, "config_flow", f"test_{os.path.basename(args.flow_config)}")
            with open(test_cfg_path, "w") as f:
                yaml.dump(config, f, default_flow_style=True)
            flow_config = test_cfg_path
            logger(f"Using test config: {test_cfg_path}", "INFO")
        else:
            flow_config = args.flow_config
    else:
        flow_config = args.flow_config

    paths = get_paths()
    n_workers = args.workers

    # ════════════════════════════════════════════════
    # Node 1: make_yaml
    # ════════════════════════════════════════════════
    if operations.get("make_yaml", False):
        logger("Step 1/4: Generating cutset YAML files...", "INFO")
        cutsets_dir = os.path.join(outdir, "cutsets")
        corr_flag = "--correlated -c" if args.correlated else ""
        cmd = f"{PYTHON} {paths['make_yaml']} {flow_config} -o {outdir} {corr_flag}"
        run_cmd(cmd)

        # Count cutset files
        mCutSets = len([f for f in os.listdir(cutsets_dir) if f.endswith(".yml")])
        logger(f"Found {mCutSets} cutset files", "INFO")
    else:
        logger("Step 1/4: make_yaml — SKIPPED", "WARNING")
        cutsets_dir = os.path.join(outdir, "cutsets")
        mCutSets = len([f for f in os.listdir(cutsets_dir) if f.endswith(".yml")]) if os.path.exists(cutsets_dir) else 0
        logger(f"Using existing {mCutSets} cutset files", "INFO")

    # ════════════════════════════════════════════════
    # Node 2: mass_fit
    # ════════════════════════════════════════════════
    if operations.get("mass_fit", False):
        logger("Step 2/4: Fitting invariant mass distributions...", "INFO")

        # Determine input: custom or auto
        custom_projs = config.get("inputs", {}).get("projs", "")
        if custom_projs:
            proj_dir = custom_projs
        else:
            proj_dir = os.path.join(outdir, "projs")

        ry_dir = os.path.join(outdir, "raw_yields")
        os.makedirs(ry_dir, exist_ok=True)

        proj_files = sorted(f for f in os.listdir(proj_dir) if f.startswith("proj_") and f.endswith(".root"))

        def run_mass_fit(proj_file):
            proj_path = os.path.join(proj_dir, proj_file)
            cmd = f"{PYTHON} {paths['mass_fit']} {flow_config} {proj_path} -o {ry_dir} -b"
            return run_cmd(cmd)

        if n_workers > 1 and len(proj_files) > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
                list(executor.map(run_mass_fit, proj_files))
        else:
            for pf in sorted(proj_files):
                run_mass_fit(pf)
    else:
        ry_dir = os.path.join(outdir, "raw_yields")
        if config.get("inputs", {}).get("raw_yields", ""):
            ry_dir = config["inputs"]["raw_yields"]
        logger("Step 2/4: mass_fit — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════
    # Node 3: cut_variation
    # ════════════════════════════════════════════════
    if operations.get("cut_variation", False):
        logger("Step 3/4: Cut variation...", "INFO")

        # Inputs: custom or auto
        custom_ry = config.get("inputs", {}).get("raw_yields", "")
        ry_path = custom_ry if custom_ry else ry_dir

        custom_eff = config.get("inputs", {}).get("effs", "")
        eff_path = custom_eff if custom_eff else os.path.join(outdir, "effs")

        cmd = f"{PYTHON} {paths['cut_variation']} {flow_config} {ry_path} {eff_path} -b"
        run_cmd(cmd)
    else:
        logger("Step 3/4: cut_variation — SKIPPED", "WARNING")

    # ════════════════════════════════════════════════
    # Node 4: data_driven_fraction
    # ════════════════════════════════════════════════
    if operations.get("data_driven_fraction", False):
        logger("Step 4/4: Data-driven fraction...", "INFO")

        custom_eff = config.get("inputs", {}).get("effs", "")
        eff_path = custom_eff if custom_eff else os.path.join(outdir, "effs")

        cutvar_path = os.path.join(outdir, "cutVar", "cutVar.root")
        if not os.path.exists(cutvar_path):
            logger(f"cutVar.root not found at {cutvar_path}, using default", "WARNING")

        cmd = f"{PYTHON} {paths['data_driven_fraction']} {cutvar_path} {eff_path} -o {outdir} -b"
        run_cmd(cmd)
    else:
        logger("Step 4/4: data_driven_fraction — SKIPPED", "WARNING")

    elapsed = time.time() - start_time
    logger(f"Analysis completed in {elapsed:.1f} seconds", "INFO")


if __name__ == "__main__":
    main()
