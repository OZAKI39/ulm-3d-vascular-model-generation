#!/usr/bin/env python3
"""Audit cached primary sources and local ColorCoded correspondence; no labeling."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import re
import json
import numpy as np
from vascular_processing.brava_branch_labels import load_exact_graph,color_source_mapping
from vascular_processing.topbrain_qc import sha256,write_json
from vascular_processing.brava_mevo_roi import write_csv


def main():
    folder=ROOT/'outputs/topbrain_brava_transfer/nn_production/type_mapping_trace'
    raw=ROOT/'vessel_model/T - Brava/swc_files/BG001.CNG.swc'
    color=raw.with_name('BG001_ColorCoded.CNG.swc')
    exact=load_exact_graph(raw);labels,mapping,version=color_source_mapping(exact,color)
    html=(folder/'BG001_official_record.html').read_text()
    mca=html.split('MIDDLE CEREBRAL ARTERY (MCA)')[1].split('MCA - Branch')[0]
    row=next(row for row in re.findall(r'<tr>(.*?)</tr>',mca,re.S) if 'Total length (mm)' in row)
    cells=re.findall(r'<td[^>]*>(.*?)</td>',row,re.S)
    lengths=[float(re.sub('<[^>]+>','',cell).strip()) for cell in cells[1:]]
    if 'LEFT</strong>' not in mca or 'RIGHT</strong>' not in mca:raise ValueError('Official table headers unavailable')
    signatures=[]
    for code in range(2,8):
        graph=exact.graph.subgraph([n for n in exact.graph if labels[n]==code])
        total=sum(np.linalg.norm(graph.nodes[a]['coords'][:3]-graph.nodes[b]['coords'][:3]) for a,b in graph.edges)
        signatures.append(dict(TYPE=code,node_count=len(graph),same_type_edge_count=graph.number_of_edges(),
            same_type_length_raw_mm=float(total),rounded_length_mm=round(float(total),2)))
    inferred={side:[r['TYPE'] for r in signatures if r['rounded_length_mm']==value] for side,value in zip(['LMCA','RMCA'],lengths)}
    downstream={'RMCA':3,'LMCA':4}
    conflict=any(inferred[side]!=[code] for side,code in downstream.items())
    # Original HTTP paths are recovered from the official record and JNLP, not
    # guessed download paths. Preserve both raw and standardized variants.
    urls={'home.php':'http://cng.gmu.edu/brava/home.php','help.php':'http://cng.gmu.edu/brava/help.php',
        'all_subjects.php':'http://cng.gmu.edu/brava/all_subjects.php?clear=1',
        'BG001_official_record.html':'http://cng.gmu.edu/brava/query_subject.php?id=1',
        'BG0002_official_record.html':'http://cng.gmu.edu/brava/query_subject.php?id=2',
        'BG001_color.jnlp':'http://cng.gmu.edu/brava/files/file_jnlp_color/BG001_color.jnlp',
        'BG001.jnlp':'http://cng.gmu.edu/brava/files/file_jnlp/BG001.jnlp',
        'BG001_official_color.swc':'http://cng.gmu.edu/brava/files/file_swc_color/BG001_ColorCoded.CNG.swc',
        'BG001_official_standardized.swc':'http://cng.gmu.edu/brava/files/file_swc/BG001.CNG.swc',
        'BG001_official_raw.swc':'http://cng.gmu.edu/brava/files/file_swc/original/BG001.CNG.swc',
        'cvapp.jar':'http://cng.gmu.edu/brava/cvapp/cvapp.jar'}
    commit='cadcae820299487a9b9ffc07c744e82065c22e8a'
    for name in ['ProcessingBravaset.py','GeneralFunctions.py','Patient.py']:
        urls[name]=f'https://raw.githubusercontent.com/rmpadmos/Coupling1DBF/{commit}/python/{name}'
    sources=[dict(file=name,url=url,sha256=sha256(folder/name)) for name,url in urls.items()]
    standard=np.loadtxt(raw);original=np.loadtxt(folder/'BG001_official_raw.swc')
    normalized=original.copy();normalized[:,3]*=-1
    report=dict(status='BRAVA_TYPE_LATERALITY_CONFLICT' if conflict else 'TYPE_MAPPING_CONSISTENT',
        subject='BG001',source_version_mapping=version,source_signatures=signatures,
        official_table=dict(source='BG001_official_record.html',html_header_lines=[523,524],html_value_lines=[544,545],
            LMCA_length_mm=lengths[0],RMCA_length_mm=lengths[1],matching_TYPEs=inferred,
            interpretation='Derived from an exact two-decimal length match, NOT an explicit TYPE legend'),
        downstream_author_code=dict(repository='https://github.com/rmpadmos/Coupling1DBF',commit=commit,
            source='GeneralFunctions.py',lines=[20,21,22,23,24,25,26],mapping=downstream,
            import_path='ProcessingBravaset.py copies SWC column 2 directly into MajorVesselID; GeneralFunctions.py names 3=R. MCA, 4=L. MCA'),
        official_color_file_byte_identical=(folder/'BG001_official_color.swc').read_bytes()==color.read_bytes(),
        official_standardized_file_byte_identical=(folder/'BG001_official_standardized.swc').read_bytes()==raw.read_bytes(),
        original_vs_standardized=dict(exact_after_y_sign_change=bool(np.array_equal(normalized,standard)),
            original_sha256=sha256(folder/'BG001_official_raw.swc'),standardized_sha256=sha256(raw),
            note='Audit only; no input coordinates changed. Local non-smoothed file is the official standardized variant.'),
        facts_not_inferences=['TYPE 3 and TYPE 4 identify MCA-associated annotated samples',
            'The numerical left/right convention is inconsistent across retrieved sources'],
        production_mapping={'LMCA':None,'RMCA':None},production_enabled=False,
        no_coordinate_based_laterality_inference=True,registration_run=False,VascularMD_run=False,
        sources=sources)
    write_json(folder/'type_mapping_audit.json',report)
    mapping_file=folder/'BG001_raw_to_color_node_mapping.csv'
    if not mapping_file.exists():write_csv(mapping_file,[dict(original_swc_id=n,color_swc_id=c,major_type=labels[n]) for n,c in mapping.items()])
    config=dict(subject='BG001',raw_source=str(raw.relative_to(ROOT)),color_source=str(color.relative_to(ROOT)),
        raw_sha256=sha256(raw),color_sha256=sha256(color),major_arteries={'LMCA':None,'RMCA':None},
        evidence_status=report['status'],evidence=str((folder/'type_mapping_audit.json').relative_to(ROOT)),
        candidate_mappings=dict(author_code=downstream,official_table_derived={k:v[0] for k,v in inferred.items()}))
    # An audit must not erase an explicitly accepted production convention.
    config_path=ROOT/'configs/brava_major_arteries.json'
    if not config_path.exists():
        write_json(config_path,config)
    print(json.dumps({k:report[k] for k in ['status','official_table','downstream_author_code','official_color_file_byte_identical','official_standardized_file_byte_identical','original_vs_standardized']},indent=2))
    return 0

if __name__=='__main__':raise SystemExit(main())
