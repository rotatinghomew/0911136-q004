"""运行领域摘要或启动本地供数服务。"""

import argparse
import json
from pathlib import Path
from http.server import ThreadingHTTPServer

from .analysis import summarize
from .datasource import SourceClient
from .lag_association import analyze_lag_association
from .reports import package_report, save_report
from .service import handler_for, load_repository, serve


def _lag_config(arguments):
    config = {}
    if arguments.lag_bin_minutes is not None:
        config["bin_minutes"] = arguments.lag_bin_minutes
    if arguments.lag_permutations is not None:
        config["permutations"] = arguments.lag_permutations
    return config


def main():
    parser = argparse.ArgumentParser(description="领域资料摘要")
    parser.add_argument("--dataset", default="fixtures/dataset.json")
    parser.add_argument("--source-url")
    parser.add_argument("--output", default="/tmp/domain-summary.json")
    parser.add_argument("--serve", action="store_true")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--study", choices=("lag-association",),
                        help="额外运行的研究；缺省只输出现行摘要")
    parser.add_argument("--lag-bin-minutes", type=int,
                        help="时滞档宽（分钟）；档宽同时是有效时滞的分辨率上限")
    parser.add_argument("--lag-permutations", type=int,
                        help="搜索校正零分布的重抽样次数")
    parser.add_argument("--lag-selection-share", type=float, default=0.5,
                        help="按过程顺序用于选型的比例；其余完整过程用于复核")
    arguments = parser.parse_args()
    schema_path = Path(__file__).with_name("schema.json")
    if arguments.serve:
        server = ThreadingHTTPServer(("127.0.0.1", arguments.port), handler_for(load_repository(arguments.dataset, schema_path)))
        try:
            print(json.dumps({"listening": server.server_address}, ensure_ascii=False), flush=True)
            server.serve_forever()
        finally:
            server.server_close()
        return
    if arguments.source_url:
        snapshot = SourceClient(arguments.source_url).read()
    else:
        with serve(load_repository(arguments.dataset, schema_path)) as base_url:
            snapshot = SourceClient(base_url).read()
    summary = summarize(snapshot.tables)
    if arguments.study == "lag-association":
        summary = {**summary, "lag_association": analyze_lag_association(
            snapshot.tables, _lag_config(arguments), arguments.lag_selection_share)}
    report = package_report(snapshot, summary)
    save_report(arguments.output, report)
    print(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
