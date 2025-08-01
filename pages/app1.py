import streamlit as st
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import os

st.set_page_config(layout="wide", page_title="DNA 돌연변이와 단백질 변화 실험")
st.title("돌연변이 실험: DNA → 단백질 변화 관찰 (단계별)")

Entrez.email = "your_email@example.com"  # 실제 이메일로 수정

# 예시 유전자 목록 (Accession, UniProt)
default_genes = [
    # --- Human genes (Accession, UniProt) ---
    ("HBB (인간 베타글로빈)",   "NM_000518", "P68871"),
    ("INS (인간 인슐린)",       "NM_000207", "P01308"),
    ("TP53 (인간 p53)",        "NM_000546", "P04637"),
    ("BRCA1 (인간 유방암 억제)", "NM_007294", "P38398"),
    ("EGFR (인간 상피세포 성장인자 수용체)", "NM_005228", "P00533"),
    # --- Other species ---
    ("lacZ (E. coli β-gal)",     "J01636", "P00722"),
    ("trpA (E. coli)",           "NP_415778", "P0A877"),
    ("COX1 (효모 mitochondria)", "NC_001224", "P00410"),
]
gene_dict = {f"{name} ({acc})": (acc, uniprot) for name, acc, uniprot in default_genes}
gene_options = list(gene_dict.keys()) + ["직접 입력"]

# 세션 상태 초기화 함수
def session_reset(start_step=0):
    keys = ["show_rna", "show_protein", "show_compare", "mut_base", "mut_pos"]
    steps = {"show_rna":1, "show_protein":2, "show_compare":3}
    for k in keys:
        if k in st.session_state:
            if start_step == 0:
                del st.session_state[k]
            else:
                # 각 단계 이상은 삭제
                if steps.get(k, 0) >= start_step:
                    del st.session_state[k]

# HTML 파일 읽기 함수
def load_html_file(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            return file.read()
    except FileNotFoundError:
        return None
    except Exception as e:
        st.error(f"HTML 파일을 읽는 중 오류 발생: {e}")
        return None

# (1) 유전자 선택 & DNA/CDS 서열 불러오기
st.markdown("#### 1️⃣ 유전자 선택/검색")
gene_choice = st.selectbox("예시 유전자/직접 입력", gene_options)
if gene_choice == "직접 입력":
    acc_input = st.text_input("GenBank Accession 번호 or RefSeq:", "")
    uniprot_id = st.text_input("(선택) UniProt ID (있을 경우 3D 구조 연동):", "")
else:
    acc_input, uniprot_id = gene_dict[gene_choice]

seq_record, cds_seq, protein_seq, cds_feature = None, None, None, None
if acc_input:
    try:
        with Entrez.efetch(db="nucleotide", id=acc_input, rettype="gb", retmode="text") as handle:
            seq_record = SeqIO.read(handle, "genbank")
        for f in seq_record.features:
            if f.type == "CDS":
                cds_feature = f
                break
        if cds_feature:
            cds_seq = cds_feature.extract(seq_record.seq)
            protein_seq = cds_feature.qualifiers.get("translation", [""])[0]
            # db_xref에서 UniProt ID 추출 (직접입력이 아니면 자동 세팅, 직접입력은 위에서)
            if gene_choice == "직접 입력" and not uniprot_id:
                for dbxref in cds_feature.qualifiers.get("db_xref", []):
                    if dbxref.startswith("UniProtKB/"):
                        uniprot_id = dbxref.split(":")[1]
                        break
    except Exception as e:
        st.error(f"NCBI에서 서열을 불러오는 중 오류 발생: {e}")

if cds_seq:
    st.success(f"CDS(코딩 DNA) 서열을 불러왔습니다! ({len(cds_seq)} bp)")
    st.code(str(cds_seq), language="text")
    st.markdown("#### 2️⃣ 원하는 위치의 염기를 바꿔보세요 (돌연변이 생성)")
    mut_col1, mut_col2 = st.columns([1,1])
    with mut_col1:
        mut_pos = st.number_input(f"돌연변이 시킬 위치 (1 ~ {len(cds_seq)})", min_value=1, max_value=len(cds_seq), value=1, key="mut_pos")
    with mut_col2:
        orig_base = str(cds_seq[st.session_state.get("mut_pos", 1)-1])
        mut_base = st.text_input(f"변경할 염기 (원래: {orig_base})", value=orig_base, key="mut_base")
        mut_base = mut_base.upper().replace(" ", "")

    is_valid = (mut_base in ["A","T","G","C"])
    mutated_seq = list(str(cds_seq))
    is_mutated = is_valid and mut_base != orig_base
    if is_valid and mut_base != orig_base:
        mutated_seq[st.session_state.get("mut_pos", 1)-1] = mut_base
        mutated_seq = "".join(mutated_seq)
        st.success(f"변이된 CDS 서열 (변이: {orig_base}{st.session_state.get('mut_pos',1)}→{mut_base})")
        st.code(mutated_seq, language="text")
    else:
        mutated_seq = str(cds_seq)

    # 단계별 세션 상태 리셋 버튼(돌연변이 수정시 아래 단계 전부 초기화) → 이곳에만 유지!
    if st.button("돌연변이 정보 수정 → 아래 과정 리셋", help="돌연변이 정보가 바뀌면 전사/번역/비교 과정을 리셋합니다."):
        session_reset(start_step=1)
        st.experimental_rerun()

    # (3) 전사하기 단계
    if st.button("전사하기 (DNA → RNA)"):
        st.session_state.show_rna = True

    if st.session_state.get("show_rna", False):
        rna_seq = str(Seq(mutated_seq).transcribe())
        st.markdown("#### ▶️ 전사 결과 (RNA)")
        st.code(rna_seq, language="text")

        # (4) 번역하기 단계
        if st.button("번역하기 (RNA → 단백질)"):
            st.session_state.show_protein = True

        if st.session_state.get("show_protein", False):
            aa_seq = str(Seq(rna_seq).translate(to_stop=False))
            st.markdown("#### ▶️ 번역 결과 (아미노산)")
            st.code(aa_seq, language="text")

            # (5) 비교 버튼
            if st.button("공식 단백질 서열과 비교"):
                st.session_state.show_compare = True

            if st.session_state.get("show_compare", False):
                def color_diff(ref, mutant):
                    out = ""
                    for a, b in zip(ref, mutant):
                        if a == b:
                            out += f"<span style='color:#228be6;font-weight:bold'>{a}</span>"
                        else:
                            out += f"<span style='color:#c92a2a;font-weight:bold'>{b}</span>"
                    if len(mutant) > len(ref):
                        out += f"<span style='color:#fa5252;font-weight:bold'>{mutant[len(ref):]}</span>"
                    return out
                st.markdown("**공식 단백질 vs 돌연변이 단백질**<br>(<span style='color:#228be6'>파란색</span>: 일치, <span style='color:#c92a2a'>빨간색</span>: 변이/불일치)", unsafe_allow_html=True)
                st.markdown(
                    f"<b>공식:</b> <span style='font-family:consolas'>{protein_seq}</span><br>"
                    f"<b>돌연변이:</b> <span style='font-family:consolas'>{color_diff(protein_seq, aa_seq)}</span>",
                    unsafe_allow_html=True
                )

                # AI 해설 프롬프트 제공
                with st.expander("AI 해설 프롬프트(복사해서 ChatGPT에 질문하세요)"):
                    prompt = f"""아래는 유전자의 CDS 염기서열에서 특정 위치가 변이된 결과와, 단백질(아미노산) 서열입니다.
- 공식 단백질: {protein_seq}
- 돌연변이 단백질: {aa_seq}
차이점(예: missense, nonsense) 및 단백질 기능 변화 예측을 중학생도 이해할 수 있게 설명해 주세요."""
                    st.code(prompt, language="text")
                    st.markdown("[ChatGPT에서 해설 보기](https://chat.openai.com/)")
                    
                # --- 공식 단백질 3D 구조(썸네일+링크) ---
                st.markdown("### 🧬 공식 단백질 3D 구조 (AlphaFold)")
                if uniprot_id:
                    st.image(
                        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4_thumbnail_medium.jpg",
                        caption="AlphaFold 3D 구조 썸네일 (야생형)"
                    )
                    st.markdown(
                        f"[AlphaFold 3D 뷰어로 이동(새 창)](https://alphafold.ebi.ac.uk/entry/{uniprot_id})"
                    )
                else:
                    st.info("UniProt ID가 없어 3D 구조 자동 연결이 불가합니다.")

                # --- 단백질 구조 예측 비교 (원본 vs 돌연변이) ---
                st.markdown("### 🧬 단백질 구조 예측 비교 (원본 vs 돌연변이)")
                with st.expander("3D 구조 비교 시각화"):
                    # HTML 파일 로드 및 표시
                    html_content = load_html_file("3dcartoon.html")
                    if html_content:
                        st.components.v1.html(html_content, height=600, scrolling=True)
                        st.write("위 3D 구조에서 원본과 돌연변이 단백질의 구조 차이를 비교해보세요.")
                    else:
                        st.warning("3dcartoon.html 파일을 찾을 수 없습니다. 파일이 같은 디렉토리에 있는지 확인해주세요.")
                        
                    # 기존 AlphaFold Colab 안내는 유지
                    st.markdown("---")
                    st.write("또는 아래 돌연변이 단백질 서열을 복사해 [AlphaFold Colab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2_mmseqs2.ipynb)에 입력하면 3D 구조를 볼 수 있습니다.")
                    st.code(aa_seq, language="text")
                    st.markdown("[AlphaFold Colab 바로가기](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2_mmseqs2.ipynb)")

else:
    st.info("위의 유전자 예시를 선택하거나 accession 번호를 입력하고, CDS 서열을 불러와주세요.")