"""Finite GUI/telemetry integration check. Close other AFG controllers first."""
import json
from pathlib import Path
import time
import tkinter as tk
import afg_gui
from afg_control import Rig, BASE

class ReadOnlyRig(Rig):
    def connect(self):
        self.commands = []
        result = super().connect()
        self.before = self.status()
        return result

    def _write(self, session, command, **kwargs):
        if not all(part.strip().endswith('?') for part in command.split(';')):
            raise RuntimeError('Hardware-setting command prohibited in read-only validation')
        self.commands.append(command)
        return super()._write(session, command, **kwargs)

    def close(self):
        self.after = self.status()
        assert not self.touched
        super().close()

if __name__ == '__main__':
    afg_gui.Rig = ReadOnlyRig
    root = tk.Tk()
    app = afg_gui.App(root)
    root.withdraw()
    began = time.monotonic()
    app.connect()
    record = {}
    def finish():
        record['before_close'] = json.loads((BASE/'observation_status.json').read_text())
        app.close()
    root.after(6500, finish)
    root.mainloop()
    record.update(wall_s=time.monotonic()-began,
                  before=app.rig.before, after=app.rig.after,
                  all_settings_unchanged=app.rig.before == app.rig.after,
                  all_query_only=all(c.strip().endswith('?') for x in app.rig.commands for c in x.split(';')),
                  hardware_control_performed=False,
                  queries=app.rig.commands)
    record['closed_status'] = json.loads((BASE/'observation_status.json').read_text())
    record['pass'] = (record['all_settings_unchanged'] and record['all_query_only']
                      and all(r['state'] == 'READ_OK' for r in record['before_close']['channels'].values())
                      and record['closed_status']['controller'] == 'CLOSED')
    (BASE/'observation_readonly_validation.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(json.dumps({k:record[k] for k in ('wall_s','all_settings_unchanged','all_query_only','pass')}))
    if not record['pass']:
        raise SystemExit(1)
