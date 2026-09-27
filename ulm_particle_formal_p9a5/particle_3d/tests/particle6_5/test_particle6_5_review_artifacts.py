import json

def test_report_and_machine_contract(report):
 d=json.loads((report/'PARTICLE6_5_VALIDATION.json').read_text());text=(report/'PARTICLE6_5_REVIEW.md').read_text()
 assert d['manual_visual_review']==d['scientific_sensitivity_review']=='PENDING_USER_REVIEW'
 assert d['near_field_regularization_v1_frozen'] and not d['particle7_started']
 assert d['raw_geometry_gap_preserved'] and d['old_p5_subnanometer_history_preserved']
 assert '## 13. 人工审核状态' in text and text.count('应该看什么')==12
 assert text.count('实际看到什么')==12 and text.count('有没有异常')==12
 assert d['effective_gap_definition']=='h_geom - h_lower'
