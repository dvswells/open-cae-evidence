"""Compatible CLI and Mission 1/2 imports; orchestration now lives in runner."""

from .adapters.calculix import check_solver_log
from .adapters.modal_dat import parse_dat
from .evidence import sha256, write_json
from .modal_report import report
from .modal_verification import TOLERANCES, convergence, evaluate, verify_modes
from .runner import main, replay, run_benchmark, run_mesh

if __name__ == "__main__":
    raise SystemExit(main())
