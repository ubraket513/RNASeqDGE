# 다음 세션 인계

> 최신 상태 (2026-09-17 후속 작업): dependency 설치와 P0 local 검증 완료,
> P1 native validator, P2 streaming count adapter/aggregation, P3 offline
> R/DESeq2 interface 구현·검증 완료. P4 alignment/dependency boundary 구현 및
> offline 검증 및 review 완료. P5 workflow 구현·검증·review 완료. P6 축소 local 비교 완료 (27m15s).
> 실제 P6 데이터는 report의 GSE80336을 사용한다 (docs/P6_REPORT_DATA.md).
> 아래는 최초 인계 기록이며, 현재 상태와
> 다음 작업은 [IMPLEMENTATION_PROGRESS.md](IMPLEMENTATION_PROGRESS.md) 및
> [LOCAL_DEVELOPMENT.md](LOCAL_DEVELOPMENT.md)를 먼저 읽는다.

2026-09-17 기준. **계획만 작성했으며 migration 구현은 아직 시작하지 않았다.**

## 사용자 결정

- C++ utilities + CPU STAR + featureCounts + 기존 R/DESeq2 방향 승인.
- HISAT2를 baseline으로 유지하며, 비교 검증 전 default 변경 없음.
- AntRepCLA의 vendored-header/Make/testing 방식을 참고.
- 이번 세션은 계획 문서 저장까지만. 다음 세션에서 구현 시작.
- 모든 분석을 C++로 한 번에 재작성하거나 R을 바로 제거하는 계획이 아님.

## 먼저 읽을 문서

1. [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md): 실행 순서·계약·검증 gate의 기준.
2. [CXX_MIGRATION.md](CXX_MIGRATION.md): 이전 조사와 dependency 선택 배경.
3. `pipeline.smk`, `merge_transcripts.py`, `metadata.py`, `run_deg_analysis.R`,
   `run_pipeline.sh`, `config.yaml.example`, `tests/test_*.sh`의 실제 코드.

## 현재 상태와 증거

- 기준 HEAD: `9f181e2310196fc22ce1159e530b320c48b4592c`.
- native source/Makefile/vendor dependencies는 아직 추가하지 않았다.
- RNASeqDGE scientific baseline, local integration, HPC benchmark는 실행하지 않았다.
- 이전 세션에서 AntRepCLA commit `806e2b0bd9d08f06c740c1e3bae80e177a3986ce`의
  36 tests/108 assertions와 golden checks를 통과했다. 별도 workers 1/2/4/8 비교도
  통과했다. 이것을 RNASeqDGE 검증 결과로 인용하지 않는다.
- `/tmp/RNASeqDGE-AntRepCLA-review`는 과거 조사 checkout이며 세션 간 존재를 보장하지 않는다.
- 마지막 Serena 재확인에서는 shell에 Node v24.21.0이 보였지만 Serena language server는
  Node PATH 오류를 반환했다. 현재도 그렇다고 단정하지 말고 재확인한다.
- `.serena/`에는 local project metadata와 memories가 있다. 그것만을 인계 근거로 삼지 않는다.
- 문서는 로컬 파일로 저장되었다. commit/push는 하지 않았다.

## 다음 세션 시작 절차

1. `git status --short`, `git rev-parse HEAD`로 변경 여부 확인. 기존 작업 보존.
2. 적용되는 AGENTS.md 확인. 라이브러리/API 사용 시 Context7 우선 조회.
3. Serena instructions를 읽고 project activation/symbol navigation 재확인.
   실패하면 원인을 보고하고 설치 여부를 추측하지 않는다.
4. 계획 P0부터 시작: tool versions/실행 가능성 조사, tiny fixture 설계,
   gene ID 및 sample metadata 계약 확인, legacy 오류 재현과 baseline 분리.
5. `docs/IMPLEMENTATION_PROGRESS.md`를 새로 만들어 각 단계의 실제 상태·명령·결과·
   blockers를 기록한다. 아직 없는 파일/명령을 구현된 것으로 보고하지 않는다.
6. P0 gate 이후 P1 native TSV/validation build를 구현한다. 기존 경로를 지우지 않는다.

## 구현 시 특히 지킬 사항

- plan에 명시된 sample/reference/design TSV schema를 먼저 확정하고 fixture로 검증.
- gene symbol annotation 실패로 gene을 버리지 않기; sample을 교집합으로 줄이지 않기.
- malformed table, integer overflow, duplicate IDs, incomplete outputs는 즉시 실패.
- DESeq2 contrast/reference를 명시하고 plot/result가 동일 비교를 사용하도록 보장.
- alignment backend 변경과 statistical-method 변경을 함께 수행하지 않기.
- compute job 안의 다운로드/패키지 설치 제거; 외부 도구 version/checksum 고정.
- Slurm resource budget, 실패 전파, atomic publish, resume invalidation 테스트.
- golden은 출력 복사만으로 정당화하지 말고 기대값·생성 환경·검토 근거 기록.

## 다음 세션에 전달할 요청문

> docs/SESSION_HANDOFF.md와 docs/IMPLEMENTATION_PLAN.md를 읽고 RNASeqDGE C++ 전환을
> P0부터 시작해줘. C++ utilities + CPU STAR + featureCounts + 기존 DESeq2가 목표이고,
> HISAT2는 검증 baseline으로 유지해. 단계별 완료 gate를 지키고 실제 진행 상황과
> 검증 결과를 docs/IMPLEMENTATION_PROGRESS.md에 기록해줘.
