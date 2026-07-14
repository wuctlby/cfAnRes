#!/usr/bin/env python3
"""
CLI entry point for the collage/merger tool.

Usage
-----
    python run.py config.yaml                  # run all tasks
    python run.py config.yaml --task NAME       # run a specific task
    python run.py config.yaml --dry-run         # validate + show plan
    python run.py config.yaml --list            # list task names
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Ensure we can import our own modules
_THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_THIS_DIR))

from schema import dry_run, load_config  # noqa: E402
from core import run_task               # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(
        description="配置驱动的拼图/PDF 合并工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python run.py config.yaml                   # 执行全部 task
  python run.py config.yaml --task my_task    # 只执行 my_task
  python run.py config.yaml --dry-run         # 干跑（校验 + 预览）
  python run.py config.yaml --list            # 列出所有 task 名称
        """.strip(),
    )
    parser.add_argument("config", help="YAML 配置文件路径")
    parser.add_argument("--task", "-t", default=None, help="只执行指定名称的 task")
    parser.add_argument("--dry-run", "-n", action="store_true", help="干跑：校验配置并预览操作，不实际生成文件")
    parser.add_argument("--list", "-l", action="store_true", help="列出所有 task 名称和描述")

    args = parser.parse_args()

    # ── load config ──
    try:
        config = load_config(args.config)
    except (FileNotFoundError, ValueError) as e:
        print(f"✗ 配置错误: {e}", file=sys.stderr)
        sys.exit(1)

    # ── list mode ──
    if args.list:
        print(f"配置文件: {args.config}")
        print(f"共 {len(config['tasks'])} 个 task:\n")
        for t in config["tasks"]:
            name = t["name"]
            ttype = t["type"]
            desc = t.get("desc", "")
            print(f"  [{ttype:11s}]  {name}")
            if desc:
                print(f"               {desc}")
        return

    # ── dry-run mode ──
    if args.dry_run:
        print(dry_run(config))
        return

    # ── run mode ──
    defaults = config.get("defaults", {})

    if args.task:
        # Find matching task
        matched = [t for t in config["tasks"] if t["name"] == args.task]
        if not matched:
            names = [t["name"] for t in config["tasks"]]
            print(f"✗ 找不到 task '{args.task}'。可用的 task: {', '.join(names)}", file=sys.stderr)
            sys.exit(1)
        tasks_to_run = matched
    else:
        tasks_to_run = config["tasks"]

    # Resolve work_dir to config's directory so relative paths in config work
    work_dir = Path(args.config).resolve().parent

    total = len(tasks_to_run)
    results: dict[str, list[Path]] = {}

    for i, task in enumerate(tasks_to_run):
        name = task["name"]
        try:
            outputs = run_task(task, defaults, work_dir=work_dir)
            results[name] = outputs
        except Exception as e:
            print(f"\n✗ Task '{name}' 失败: {e}", file=sys.stderr)
            if i < total - 1:
                print("  继续执行后续 task...\n")
            else:
                sys.exit(1)

    # ── summary ──
    print(f"\n{'=' * 60}")
    print(f"完成! 共执行 {total} 个 task")
    for name, outputs in results.items():
        print(f"  [{name}] → {len(outputs)} 个输出文件")
        for op in outputs[:5]:  # show first 5
            print(f"           {op}")
        if len(outputs) > 5:
            print(f"           ... 还有 {len(outputs) - 5} 个文件")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
