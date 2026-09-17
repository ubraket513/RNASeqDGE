# RNASeqDGE C++ 전환 실행 계획

작성일: 2026-09-17

> 실행 진행 상황은 [IMPLEMENTATION_PROGRESS.md](IMPLEMENTATION_PROGRESS.md)가
> 최신 기준이다. 아래의 "미시작" 표기는 계획 작성 당시 상태이다.

상태: **계획 수립 완료, 구현 미시작**. 사용자는 이 문서를 먼저 저장하고
다른 세션에서 구현을 시작하기로 했다. 이 문서의 명령과 파일 구조는 별도
표시가 없는 한 구현 목표이며, 현재 존재하는 기능으로 해석하면 안 된다.

기준 Git HEAD: `9f181e2310196fc22ce1159e530b320c48b4592c`.
다음 세션에서 HEAD와 작업 트리를 다시 확인하고 사용자 변경을 보존한다.

## 1. 승인된 방향과 범위

첫 제품 목표는 **C++ utilities + CPU STAR + featureCounts + 기존 DESeq2**다.
STAR는 우선 평가 대상이며, 검증 전에는 HISAT2 기본값을 변경하지 않는다.
기존 HISAT2 경로를 비교 기준과 복구 경로로 유지한다.

- 자체 유틸리티: C++20, g++, GNU Make, vendored header libraries.
- 통계: R/DESeq2/apeglm 유지. 실행 시 패키지를 설치하지 않는 고정 환경.
- 오케스트레이션: 로컬 Make와 Slurm 제출 계층. 단계적으로 Snakemake 대체.
- 정렬·counting: 검증된 native 실행 파일을 subprocess로 호출.
- metadata/reference 준비와 실제 분석 실행을 분리. 준비된 입력으로 분석은
  네트워크 연결 없이 실행 가능해야 한다.
- AntRepCLA의 vendoring 및 golden-output 검증 방식을 참고한다.
- 다른 파이프라인과 공유할 코드는 TSV, 입력 검증, 프로세스 실행, provenance
  등 실제 중복이 확인된 기능부터 추출한다. 범용 workflow engine을 새로 만들지 않는다.

이번 범위에서 제외: native DESeq2 재구현, GPU/Parabricks 개선, 자동 통계모형
선택, 모든 포맷의 직접 구현, 기존 연구 결과의 무조건적인 재현 주장.
Salmon/fastp/minimap2는 후속 평가 대상으로 두고 첫 구현 의존성에 넣지 않는다.

## 2. 현재 구현과 확인할 위험

| 파일 | 현재 역할 | 이전 시 처리 |
| --- | --- | --- |
| `pipeline.smk` | 다운로드, 정렬, counting, merge, DEG DAG | 검증용 legacy 경로로 보존 후 교체 |
| `run_pipeline.sh` | config와 job script 생성, 내부 sbatch 호출 | 실행별 작업 디렉토리와 단일 제출 계층으로 교체 |
| `merge_transcripts.py` | BioMart annotation과 count inner join | reference gene ID 기반 native merge로 교체 |
| `metadata.py` | SRR→GSM 조회와 cache | staging 단계로 이동, 명시적 mapping 입력 사용 |
| `run_deg_analysis.R` | GEO metadata, DESeq2, plotting | 오프라인 count/sample/annotation 인터페이스로 수정 |
| `requirements.txt` | unpinned Python dependencies | legacy 유지 기간 후 새 경로에서 제거 |
| `tests/test_*.sh` | cluster smoke tests | 파일 존재 검사를 넘어 값·실패·resume 검증 추가 |

기존 코드에서 확인된 사항:

- 매 job마다 임시 Python 환경에 패키지를 설치하며 Snakemake는 requirements에 없다.
- gene_id를 gene symbol로 해석하는 경로와 BioMart 실패 시 gene 누락이 있다.
- GEO와 count의 교집합만 사용하여 sample이 제외될 수 있다.
- condition은 title에서 추론하고 reference는 알파벳순이다.
- 결과와 shrinkage가 동일 contrast를 가리키는지 명시적 검증이 필요하다.
- HISAT2와 동시 실행되는 samtools sort에 각각 전체 thread budget을 준다.
- launcher가 기존 config를 덮어쓴다. README의 외부 sbatch와 내부 sbatch가 겹친다.
- counting summary를 삭제하며 strandedness를 명시하지 않는다.
- `rule all`의 기본 target 지정, awk의 tab escaping, condition 정규식 escaping,
  fasterq-dump의 read 제한 옵션은 fixture와 실제 고정 버전에서 확인해야 한다.
  아직 실행 검증하지 않은 항목을 확정 버그로 취급하지 않는다.

기존 pipeline이 실패하면 실패 로그를 남긴다. 잘못된 결과를 golden으로 고정하지
않는다. 원본 실행 결과와 검토 후 수정한 scientific baseline을 별도로 보관한다.

## 3. 구조와 의존성 결정

예정 구조:

```text
src/                  C++ CLI, tables, manifests, counts, process, provenance
include/rnaseq/        프로젝트 내부 인터페이스
third_party/          고정된 작은 header 라이브러리와 notices
vendor/               필요한 upstream 도구 source archives와 manifest
tools/                도구 준비, local/Slurm 실행, 검증 scripts
config/               sample/design/reference/toolchain 예제
tests/unit/           C++ 단위 테스트
tests/fixtures/       작은 합성 reads/reference/counts/metadata
tests/golden/         출처가 기록된 검증 결과
tests/integration/    실행, 실패, resume 테스트
docs/                 계획, 사용법, 결정, benchmark 보고서
build/                빌드 결과; Git 제외
runs/<run_id>/         실행별 immutable 입력 snapshot, 결과, 로그
```

### 자체 C++ 빌드

- csv-parser 채택. AntRepCLA에서 확인한 5.3.0은 후보 pin이며 upstream 출처,
  checksum, compiler 호환성을 확인한 후 확정한다.
- doctest 채택. 관측 header의 2.5.0 표기만으로 release provenance를 판단하지 않는다.
- 기본 container는 STL. unordered_dense는 profiling 근거가 있을 때 추가한다.
- 처음에는 별도 logging/CLI/plotting library를 추가하지 않는다.
- OpenMP는 첫 merge 기능의 필수 의존성이 아니다. 측정된 병렬 작업에만 도입한다.
- 자체 유틸리티 기준 toolchain은 우선 GCC 12 이상/C++20을 목표로 하며,
  vendor 호환성을 실제 빌드로 확인하고 최소 지원 버전을 기록한다.
- `.d` header dependencies, 사용자 지정 CXX/CPPFLAGS/CXXFLAGS/LDFLAGS를 지원한다.
  기본 `-march=native`, host 전체 core 강제 사용, 빌드 중 다운로드는 금지한다.
- native 기본 `make check`는 R, Slurm, 인터넷 없이 실행 가능하게 한다.

### 외부 분석 도구

STAR, HISAT2, featureCounts, samtools는 각각 독립 실행 파일로 관리한다.
source archive/commit, SHA-256, license/notices, build flags, 실제 version 출력을
toolchain manifest에 기록한다. 외부 library들도 transitive dependency 목록에 포함한다.

`make`는 사용자 진입점이지 모든 upstream 도구가 GCC/Make만으로 빌드된다는
보장은 아니다. native 빌드와 R 환경 준비의 요구사항을 분리해서 문서화한다.
도구 준비는 명시적 별도 단계이며, compute job에서 fetch/install하지 않는다.

HISAT2 baseline은 현재 설정의 2.2.1로 시작한다. 해당 버전의 alignment wrapper는
Perl, build wrapper는 Python이므로 legacy 경로는 interpreter-free가 아니다.
첫 단계에서 wrappers를 재작성하지 않는다. CPU STAR 신규 경로로 interpreter
의존성을 줄이되 R backend는 유지한다. upstream source와 라이선스를 보존하고
도구 내부 코드를 자체 유틸리티에 직접 합치지 않는다.

R은 선택된 R/Bioconductor와 모든 패키지 버전을 고정한다. 실제 preflight는
library 로드, tiny DESeq2/apeglm 실행, plotting device까지 검사한다. 환경 구성
기록과 sessionInfo를 보존한다. 컨테이너는 지원 환경에서의 선택적 배포 방식이다.

## 4. 데이터 규약 v1

### 공통 표 규칙

UTF-8 TSV, 첫 행 header, LF 출력/CRLF 입력 허용. BOM은 파일 시작에서만 처리한다.
column 이름은 대소문자 구분, 중복 금지. tab delimiter를 지정하고 추론하지 않는다.
v1 manifest의 셀 내부 tab/CR/LF/NUL은 금지하며 CSV-style quoted field는 지원한다.
주석 행은 manifest에서 허용하지 않는다. featureCounts의 주석 header는 전용
adapter에서 처리한다. ragged row는 즉시 오류로 처리한다.

CSV parser에는 `VariableColumnPolicy::THROW`를 사용하고, headerless legacy
counts는 별도의 field-count 검증을 한다. 필수 셀, 중복 ID, 숫자 overflow,
존재하지 않는 파일은 filename/record/column 정보를 포함하여 실패시킨다.

### 입력 파일

| 파일 | 필수 필드 / 규칙 |
| --- | --- |
| `samples.tsv` | `sample_id`, `condition`; optional covariates는 설계에서 선언 |
| `runs.tsv` | `run_id`, `sample_id`, `fastq_1`, `fastq_2`, `layout`, `strandedness` |
| `references.tsv` | `role`, `path`, `sha256`, `source`, `release`; genome/annotation 각 1개 |
| `contrasts.tsv` | `contrast_id`, `factor`, `numerator`, `denominator` |
| `analysis.tsv` | `key`, `value`; version, design_terms, alpha, filter, shrinkage 등 고정 설정 |

sample/run/contrast ID는 `[A-Za-z0-9][A-Za-z0-9_.-]*`, `.`과 `..` 단독은 금지한다.
sample은 실제 분석 단위이며 SRR/GSM 이름을 강제하지 않는다. 여러 run이 같은
sample을 참조할 수 있으나 sample metadata 자체는 unique여야 한다.

`layout`은 `single`/`paired`, `strandedness`는 `unstranded`/`forward`/`reverse`다.
single의 fastq_2는 빈 값, paired는 두 경로가 필수다. 경로는 manifest 디렉토리를
기준으로 해석하고 실행 snapshot에는 resolved path와 checksum을 기록한다.
추론한 library 속성은 후보로 제시할 수 있지만 조용히 확정하지 않는다.

기술 run을 sample로 합치는 것은 `runs.tsv`의 명시적 매핑에 의해서만 수행한다.
초기 구현은 같은 sample의 layout/strandedness가 일치해야 한다. 불일치하면
실패하고 별도 방법 검토를 요구한다. biological replicate를 합치지 않는다.

### Counts와 annotation

- `counts.tsv`: `gene_id` 다음으로 `samples.tsv` 순서의 sample columns.
- gene ID는 reference의 원문을 보존한다. symbol 변환이나 version suffix 제거를
  자동 적용하지 않는다. gene 순서는 canonical reference 목록의 bytewise 정렬.
- raw featureCounts와 `.summary`를 보존하고 전용 adapter가 canonical TSV를 생성한다.
- 같은 reference gene set이 모든 run에 존재해야 한다. 누락은 0으로 보정하지 않는다.
- annotation은 `gene_id`, `gene_symbol`, `biotype`, `chromosome`의 별도 표다.
  optional annotation 누락은 허용하지만 counts를 제거하거나 행을 증식시키지 않는다.
- 합계는 checked uint64 연산. R 전달 전 DESeq2의 지원 정수 범위와 변환을 검증하고
  범위를 벗어난 값은 silently coerce하지 않는다.
- 기술 run별 counting 후 sample별 합산한다. 기존 duplicate 처리 없는 counting
  정책을 기준으로 하며, deduplication 도입 시 이 집계 정책도 재검토한다.

### 통계 경계

R은 count/sample/annotation 파일과 명시적 analysis/contrast 설정만 읽는다.
GEO/BioMart 요청, title parsing, sample 교집합 선택을 제거한다.

v1 design은 선언된 변수의 additive model만 지원한다. arbitrary R expression을
실행하지 않는다. categorical level, numeric covariate, missing value, rank deficiency,
residual degrees of freedom과 contrast 존재 여부를 검사한다. paired design은
명시적 donor 변수를 통해 표현하고 식별 가능성을 검증한다.

초기 parity 설정: zero-total gene 제외, DESeq2 기본 normalization/dispersion,
alpha 0.05, 명시적 condition 비교, apeglm shrinkage. 사용한 coefficient와
unshrunken 결과의 방향을 일치시킨다. 복잡한 contrast의 shrinkage는 v1에서
지원하지 않고 명시적으로 오류 처리한다.

FDR cutoff는 표·로그·volcano에 일관되게 사용한다. top-ranked heatmap과 significant
gene heatmap을 구분하고 DEG 0개인 경우를 정상 처리한다. 작은 fixture에서 VST나
plot이 불가능하면 그 이유를 명시하며 임의의 통계값으로 대체하지 않는다.

## 5. 실행·재시작·HPC 설계

흐름: input validation → reference/index → run alignment → run counting →
sample merge → DESeq2 → report/provenance.

- CPU STAR alignment와 indexing은 신규 backend다. 기존 Parabricks branch와 분리한다.
  index 설정은 reference/GTF/STAR version/read-length 정책을 포함하여 고정한다.
- 초기에는 STAR의 unsorted BAM 이후 별도 samtools sort 단계로 resource accounting을
  단순화한다. direct sorted BAM 최적화는 측정 이후 수행한다.
- paired fragment counting, strandedness, MAPQ, multimapping, overlap 정책을 명시한다.
  aligner 변경과 counting 설정 변경은 같은 benchmark에서 동시에 수행하지 않는다.
- 전역 root config를 덮어쓰지 않고 `runs/<run_id>/`별로 설정을 snapshot한다.
- subprocess는 argv를 분리하여 전달한다. 사용자 경로/값을 shell 문자열로 결합하지 않는다.
  종료 코드, stderr, signal forwarding, child 종료를 처리한다.
- thread 수는 명시적 설정을 우선하며 Slurm allocation을 초과하면 실패시킨다.
  기본 local worker는 1. parser/OpenMP/BLAS/도구별 auxiliary thread도 budget에 포함한다.
- Slurm은 reference 준비 job → alignment/counting arrays → merge/DE job으로 나눈다.
  `afterok` 의존성과 array concurrency limit을 사용하고 job ID를 기록한다.
- local Make의 `-j`만으로 heterogeneous task memory를 관리하지 않는다. 최초 구현은
  resource-heavy 단계의 동시성을 보수적으로 제한한다.
- scratch는 실행·task별 전용 경로를 사용하고 성공적으로 publish된 후에만 정리한다.
  실패 로그와 재현 정보는 유지한다.
- output은 같은 filesystem의 임시 경로에서 검증 후 atomic rename한다. 여러 파일은
  전체 검증 후 completion manifest를 마지막에 게시한다. 파일 존재만으로 skip하지 않는다.
- stage identity에는 input hashes, tool versions, args, reference/config hashes가 들어간다.
  수정된 upstream stage의 downstream 결과는 invalidate한다. mtime만으로 판단하지 않는다.
- 동일 run directory의 동시 실행을 막는다. 중단·실패 시 completion 표식을 남기지 않는다.
- reference index cache는 버전별 key와 lock으로 보호하고 완성 후 publish한다.
- `plan`/dry-run은 실행할 단계와 자원을 출력하며 다운로드나 job 제출을 하지 않는다.

## 6. 단계별 구현 순서와 완료 조건

모든 단계는 **미시작**이다. 각 단계 종료 시 progress 문서에 변경 파일, 실행 명령,
실제 결과, 미해결 항목을 기록한다. 시간 추정 대신 아래 gate로 완료를 판단한다.

| 단계 | 작업 / 산출물 | 완료 gate |
| --- | --- | --- |
| P0 | toolchain 점검, legacy 오류 재현, 작은 fixture와 baseline 출처 기록 | 입력·gene/sample mapping 검토, 원본/수정 baseline 분리, 재실행 가능 |
| P1 | native Make build, vendored csv-parser/doctest, TSV validation CLI | 오프라인 build/check 통과, malformed/duplicate/overflow 테스트 |
| P2 | counts adapter/merge, technical-run aggregation, annotation 분리 | 수작업 oracle 및 검토된 baseline과 integer matrix 완전 일치 |
| P3 | R 오프라인 인터페이스, 명시적 design/contrast, 환경 고정 | sample 누락 없이 동일 입력 통계 parity, NA/방향/plot 의미 확인 |
| P4 | tool manifest/build staging, CPU STAR와 HISAT2 backend | tiny SE/PE fixture 정렬·counting, 실패 전파, thread/resource 제한 검증 |
| P5 | local Make workflow와 Slurm arrays, resume/provenance | 취소/재시작/부분 output/변경 입력/동시 실행 테스트 통과 |
| P6 | real subset 및 대표 full data에서 HISAT2↔STAR 평가 | 성능·과학적 차이 보고서, 기본 backend 결정 근거 기록 |
| P7 | 문서/CI/release 정리와 검증된 경로 전환 | 새 checkout에서 재현, legacy 복구 경로, 사용자 사용법 완성 |

P0→P1→P2→P3 순서로 scientific input 경계를 먼저 확정한다. P4는 P0 이후 독립적으로
준비할 수 있으나 P5 통합은 P2–P4 gate 통과 후다. P6 이전에 default를 바꾸지 않는다.
legacy Python/Snakemake 파일의 삭제는 P7 gate 통과 후 별도 검토한다.

## 7. 검증과 benchmark 기준

### 구현 이전의 동등성

- 작은 합성 FASTA/GTF/SE·PE reads에는 알려진 exon/junction/strand/ambiguous cases를 넣는다.
- counts oracle은 수작업으로 검토 가능한 작은 정수 표를 사용한다.
- gene ID, count, sample 순서와 NA/filter masks는 정확히 비교한다.
- 동일 pinned R 환경·동일 입력의 통계 비교 초기 허용 오차는
  `abs(a-b) <= 1e-10 + 1e-7*abs(reference)`로 설정한다. 극소 p-value에는
  log-scale 비교를 추가하고 0/NA 처리 규칙을 기록한다. 임계값 근처 DEG membership
  변경은 별도 보고한다. 실패를 숨기려고 tolerance를 넓히지 않는다.
- PNG byte parity를 요구하지 않는다. plot 입력 데이터, labels, contrast, 출력 생성과
  렌더링을 검증한다. raw/shrunken LFC를 구별한다.
- workers 1/2/4와 지원되는 다른 수를 **실제로 전달**해서 결과를 비교한다.
- 표준 test에는 download, 전체 인간 reference, cluster 할당을 요구하지 않는다.
  native unit, R integration, tool integration, real-data HPC test를 구분한다.
- sanitizer 빌드로 parser/merge/process 경계의 메모리·undefined behavior를 점검한다.

### STAR와 HISAT2의 방법 비교

같은 reads, annotation, counting 정책, sample design, R 환경을 사용한다.
정렬 방법이 다르므로 BAM/count/DEG의 완전 일치를 합격 기준으로 요구하지 않는다.

측정: stage/end-to-end wall time, CPU time, peak RSS, scratch/disk 사용량,
mapping/unique/multimapping rates, junction 및 assignment summary, gene별 count
차이, LFC 방향/크기, 유의 gene 집합 차이. 비교 가능한 budget으로 최소 3회 실행해
분산을 기록하고 index 구축 비용과 재사용 실행, cold/warm cache 조건을 구분한다.

실제 데이터는 첫 N개 SRR이 아닌 condition/replicate를 보존한 명시적 subset을 쓴다.
전체 실행 전 두 조건과 필요한 covariate가 존재하고 model이 식별 가능한지 확인한다.
속도만으로 STAR를 승격하지 않는다. synthetic expected behavior 통과와 실제
불일치 원인 검토가 필수이며 설명되지 않은 systematic bias가 있으면 HISAT2를 유지한다.

## 8. 후속 평가와 미확정 항목

- Salmon: 별도 transcript quantification route. tximport 등 적절한 import와
  effective-length 처리를 포함해야 하며 TPM을 raw counts처럼 전달하지 않는다.
- fastp: QC/adapter 필요성과 read 보존 정책을 먼저 결정한다. 무조건 trimming하지 않는다.
- HTSlib: 직접 BAM 처리가 필요한 기능이 생길 때만 링크한다.
- SeqAn3: 자체 서열 알고리즘 요구가 생길 때 검토한다.
- minimap2: long-read 및 새로운 short-read splice backend 평가 후보다.

다음 세션에서 확인할 사항과 기본 진행 방침:

| 미확정 정보 | 확인 시점 / 방침 |
| --- | --- |
| cluster OS/CPU, Slurm partition/QoS, scratch와 memory limit | P0 조사; 모르면 native local 단계는 계속 진행 |
| toolchain 실제 설치와 정확한 release pins | P0/P4; 현재 module 이름을 사용 가능성으로 간주하지 않음 |
| sample/run mapping과 strand, covariates | P0 scientific baseline 전에 실제 metadata 검증 |
| full dataset 접근과 분석 resource budget | P6 전에 확정; synthetic/local 테스트로 대체 완료 처리 금지 |
| R 배포 환경 방식 | P3 전 preflight 가능한 고정 환경 선택 |

질문이 필요한 경우 위 정보에 의존하는 단계만 보류한다. 추측으로 biological
metadata를 채우거나 cluster 자원을 무제한 사용하지 않는다.

## 9. 출처와 문서 관계

- [CXX_MIGRATION.md](CXX_MIGRATION.md): AntRepCLA 조사, 관측 버전과 테스트 증거.
- [SESSION_HANDOFF.md](SESSION_HANDOFF.md): 다음 세션의 시작점과 현재 상태.
- [AntRepCLA 고정 commit](https://github.com/ubraket513/AntRepCLA/tree/806e2b0bd9d08f06c740c1e3bae80e177a3986ce)
- [STAR source/build](https://github.com/alexdobin/STAR)
- [nf-core RNA-seq usage](https://github.com/nf-core/rnaseq/blob/master/docs/usage.md)
- [HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/)
- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [Salmon build requirements](https://salmon.readthedocs.io/en/latest/building.html)

이 문서는 구현 순서와 범위의 기준이다. 기존 조사 문서와 충돌하면 이 문서의
최신 결정이 우선한다. upstream API/버전은 구현 시 Context7로 먼저 조회하고,
없거나 부정확하면 고정 upstream source와 공식 문서로 확인한다.
