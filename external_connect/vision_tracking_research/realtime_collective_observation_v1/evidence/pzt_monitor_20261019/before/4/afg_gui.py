"""English desktop UI; hardware operations run on one background worker."""
import json
import logging
from pathlib import Path
import queue
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict

from afg_control import BASE, ChannelSettings, Rig, load_config

MODES = {"Continuous Sine": "continuous", "Sweep": "sweep", "N-Cycle Burst": "burst", "External Gate": "gated"}
SOURCES = {"Manual / SW (MAN)": "MAN", "Internal (TIM)": "TIM", "External (EXT)": "EXT"}


def configure_dpi():
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        except (AttributeError, OSError):
            pass


class App:
    def __init__(self, root):
        self.root = root
        self.config = load_config()
        self.rig = None
        self.executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix="AFG-USB")
        self.events = queue.Queue()
        self.busy = False
        self.closing = False
        self.buttons = []
        self.fields = {}
        self.states = {}
        self.last_readback = {}
        root.title("AFG1062 - Four-Channel PZT Control")
        root.geometry("1240x850")
        root.minsize(1000, 740)
        style = ttk.Style()
        style.theme_use("clam")
        style.configure("TLabel", font=("Segoe UI", 10))
        style.configure("TButton", font=("Segoe UI", 10), padding=5)
        style.configure("Title.TLabel", font=("Segoe UI", 17, "bold"))
        style.configure("Stop.TButton", foreground="#a51d25", font=("Segoe UI", 11, "bold"))
        fixed = ttk.Frame(root, padding=(16, 6))
        fixed.pack(fill="x")
        ttk.Label(fixed, text="AFG1062 Output Control").pack(side="left")
        ttk.Button(fixed, text="All Outputs OFF", style="Stop.TButton", command=self.stop).pack(side="right")
        # Scrollbars keep every control reachable on smaller / high-DPI screens.
        viewport = ttk.Frame(root)
        viewport.pack(fill="both", expand=True)
        viewport.rowconfigure(0, weight=1)
        viewport.columnconfigure(0, weight=1)
        canvas = tk.Canvas(viewport, highlightthickness=0)
        canvas.grid(row=0, column=0, sticky="nsew")
        vertical = ttk.Scrollbar(viewport, orient="vertical", command=canvas.yview)
        vertical.grid(row=0, column=1, sticky="ns")
        horizontal = ttk.Scrollbar(viewport, orient="horizontal", command=canvas.xview)
        horizontal.grid(row=1, column=0, sticky="ew")
        canvas.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
        outer = ttk.Frame(canvas, padding=16)
        item = canvas.create_window((0, 0), window=outer, anchor="nw")
        def resize(_=None):
            canvas.itemconfigure(item, width=max(outer.winfo_reqwidth(), canvas.winfo_width()),
                                 height=max(outer.winfo_reqheight(), canvas.winfo_height()))
            canvas.configure(scrollregion=canvas.bbox("all"))
        outer.bind("<Configure>", resize)
        canvas.bind("<Configure>", resize)
        ttk.Label(outer, text="Two Generators - Four Bridge Amplifiers - Four PZTs", style="Title.TLabel").pack(anchor="w")
        ttk.Label(outer, text="Enter generator output / amplifier input in Vpp (50 Ω). Bridge PZT voltage is approximately 14x higher; measure the actual voltage.",
                  foreground="#805000").pack(anchor="w", pady=(4, 10))
        bar = ttk.Frame(outer)
        bar.pack(fill="x")
        self.button(bar, "Connect Devices (Read Only)", self.connect).pack(side="left", padx=(0, 6))
        self.button(bar, "Read Instrument Status", self.refresh).pack(side="left", padx=6)
        self.button(bar, "Apply All (Outputs OFF)", self.apply_all).pack(side="left", padx=6)
        self.connection = tk.StringVar(value="Disconnected. Outputs are never enabled automatically.")
        ttk.Label(outer, textvariable=self.connection, wraplength=1150).pack(anchor="w", pady=8)
        cards = ttk.Frame(outer)
        cards.pack(fill="both", expand=True)
        for index, (name, ch) in enumerate(self.config["channels"].items()):
            cards.columnconfigure(index, weight=1, uniform="channel")
            card = ttk.LabelFrame(cards, text=f"{name}  /  {ch['device']} · CH{ch['channel']}", padding=10)
            card.grid(row=0, column=index, sticky="nsew", padx=4)
            card.columnconfigure(1, weight=1)
            self.states[name] = tk.StringVar(value="Output: unknown")
            ttk.Label(card, textvariable=self.states[name], foreground="#175075", wraplength=250).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 8))
            ttk.Label(card, text=f"Input limit: {ch['max_input_vpp']:g} Vpp\nFrequency limit: {ch['max_frequency_hz']/1e6:g} MHz",
                      foreground="#805000").grid(row=1, column=0, columnspan=2, sticky="w", pady=(0, 8))
            self.fields[name] = {}
            entries = [
                ("frequency", "Frequency / MHz", "2.5", None),
                ("amplitude", "Input / Vpp", "0.1", None),
                ("phase", "Phase / °", "0", None),
                ("mode", "Mode", "Continuous Sine", list(MODES)),
                ("start", "Sweep Start / MHz", "2", None),
                ("stop", "Sweep Stop / MHz", "3", None),
                ("time", "Sweep Time / ms", "100", None),
                ("spacing", "Sweep Spacing", "LIN", ["LIN", "LOG"]),
                ("cycles", "Burst Cycles", "10", None),
                ("source", "Trigger Source", "Manual / SW (MAN)", list(SOURCES)),
            ]
            for row, (key, title, initial, values) in enumerate(entries, start=2):
                ttk.Label(card, text=title).grid(row=row, column=0, sticky="w", pady=5)
                var = tk.StringVar(value=initial)
                self.fields[name][key] = var
                if values:
                    widget = ttk.Combobox(card, textvariable=var, values=values, state="readonly", width=14)
                else:
                    widget = ttk.Entry(card, textvariable=var, width=14)
                widget.grid(row=row, column=1, sticky="ew", pady=5, padx=(6, 0))
                if key in ("frequency", "amplitude", "phase"):
                    widget.bind("<Return>", lambda event, n=name: self.live_one(n))
                    widget.bind("<KP_Enter>", lambda event, n=name: self.live_one(n))
            rowbar = ttk.Frame(card)
            rowbar.grid(row=12, column=0, columnspan=2, sticky="ew", pady=(12, 5))
            self.button(rowbar, "Apply (OFF)", lambda n=name: self.apply_one(n)).pack(side="left", expand=True, fill="x")
            self.button(rowbar, "Apply Live", lambda n=name: self.live_one(n)).pack(side="left", padx=(4, 0), expand=True, fill="x")
            switches = ttk.Frame(card)
            switches.grid(row=13, column=0, columnspan=2, sticky="ew", pady=(0, 5))
            self.button(switches, "ON", lambda n=name: self.output(n, True)).pack(side="left", padx=(0, 4), expand=True, fill="x")
            self.button(switches, "OFF", lambda n=name: self.output(n, False)).pack(side="left", expand=True, fill="x")
            ttk.Label(card, text="Apply Live / Enter: frequency, Vpp, phase.\nContinuous Sine; keeps ON/OFF state.\nFirst setup or mode change: Apply (OFF).",
                      wraplength=250, foreground="#555555").grid(row=14, column=0, columnspan=2, sticky="w", pady=8)
        cards.rowconfigure(0, weight=1)
        lower = ttk.Frame(outer)
        lower.pack(fill="x", pady=10)
        for alias in ("A", "B"):
            self.button(lower, f"{alias}: Align CH1/CH2", lambda a=alias: self.align(a)).pack(side="left", padx=4)
            self.button(lower, f"{alias}: Trigger MAN", lambda a=alias: self.trigger(a)).pack(side="left", padx=4)
        self.button(lower, "Save Settings", self.save_profile).pack(side="right", padx=4)
        self.button(lower, "Load Settings", self.load_profile).pack(side="right", padx=4)
        ttk.Label(outer, text="A shared 10 MHz clock requires cabling, front-panel setup and phase calibration. USB commands do not ensure simultaneous output.\n"
                  "All Outputs OFF may fail if USB communication is lost. Normal exit attempts to turn off devices controlled in this session.",
                  foreground="#555555").pack(anchor="w")
        self.log = tk.Text(outer, height=6, wrap="word", state="disabled", font=("Segoe UI", 9))
        self.log.pack(fill="x", pady=(8, 0))
        root.update_idletasks()
        root.geometry(f"{min(outer.winfo_reqwidth()+32, root.winfo_screenwidth()-80)}x"
                      f"{min(outer.winfo_reqheight()+32, root.winfo_screenheight()-100)}+20+20")
        root.protocol("WM_DELETE_WINDOW", self.close)
        root.after(20, self.poll)

    def button(self, parent, label, action):
        button = ttk.Button(parent, text=label, command=action)
        self.buttons.append(button)
        return button

    def note(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def require_rig(self):
        if self.rig is None:
            raise RuntimeError("No devices connected. Connect the instruments first.")
        return self.rig

    def submit(self, title, function, callback=None):
        if self.busy or self.closing:
            return
        self.busy = True
        for button in self.buttons:
            button.state(["disabled"])
        self.note(title + "…")
        future = self.executor.submit(function)
        future.add_done_callback(lambda f: self.events.put((title, f, callback)))

    def poll(self):
        try:
            while True:
                title, future, callback = self.events.get_nowait()
                self.busy = False
                for button in self.buttons:
                    button.state(["!disabled"])
                try:
                    result = future.result()
                    if callback:
                        callback(result)
                    self.note(title + ": complete")
                except Exception as e:
                    logging.exception(title)
                    self.note(title + ": failed - " + str(e))
                    for state in self.states.values():
                        state.set("Output: read status again")
                    self.last_readback.clear()
                    if not self.closing:
                        messagebox.showerror(title, str(e))
                if self.closing:
                    self.finish_close()
                    return
        except queue.Empty:
            pass
        self.root.after(20, self.poll)

    def connect(self):
        if self.rig is not None:
            self.refresh()
            return
        def action():
            rig = Rig(self.config)
            rig.connect()
            return rig
        def done(rig):
            self.rig = rig
            self.connection.set("Connected: " + "    |    ".join(f"{a} = {d['serial']}  Firmware {d['firmware']}" for a, d in rig.identities.items()))
            self.note("Connection read instrument identities only. No channel settings or outputs were changed.")
        self.submit("Connect Devices", action, done)

    def collect(self, name):
        f = {k: v.get().strip() for k, v in self.fields[name].items()}
        mode = MODES[f["mode"]]
        result = ChannelSettings(
            frequency_hz=float(f["frequency"]) * 1e6,
            amplitude_vpp=float(f["amplitude"]), phase_deg=float(f["phase"]), mode=mode,
            sweep_start_hz=float(f["start"]) * 1e6, sweep_stop_hz=float(f["stop"]) * 1e6,
            sweep_time_s=float(f["time"]) / 1000, sweep_spacing=f["spacing"],
            trigger_source="EXT" if mode == "gated" else SOURCES[f["source"]], burst_cycles=int(f["cycles"]),
        )
        result.validate(self.config["channels"][name])
        return result

    def show_states(self, result):
        for name, value in result.items():
            if name in self.states:
                previous = self.last_readback.setdefault(name, {})
                previous.update(value)
                text = f"Output: {'ON' if value['output'] else 'OFF'}"
                if "frequency_hz" in previous:
                    text += (f"\nLast read: {previous['frequency_hz']/1e6:g} MHz / {previous['amplitude_vpp']:g} Vpp"
                             f"\nPhase {previous['phase_deg']:.2f}°")
                self.states[name].set(text)

    def refresh(self):
        self.submit("Read Instrument Status", lambda: self.require_rig().status(), self.show_states)

    def apply_names(self, names):
        try:
            requested = {name: self.collect(name) for name in names}
        except (ValueError, KeyError) as e:
            messagebox.showerror("Invalid Settings (Nothing Sent)", str(e))
            return
        self.submit("Apply " + ", ".join(names), lambda: self.require_rig().apply(requested), self.show_states)

    def apply_one(self, name):
        self.apply_names([name])

    def apply_all(self):
        self.apply_names(list(self.fields))

    def live_one(self, name):
        if self.busy or self.closing:
            if not self.closing:
                self.note("Busy: wait for the current operation, then Apply Live again.")
            return "break"
        try:
            f = self.fields[name]
            if MODES[f["mode"].get()] != "continuous":
                raise ValueError("Apply Live supports Continuous Sine only. Use Apply (OFF) for mode changes.")
            settings = ChannelSettings(frequency_hz=float(f["frequency"].get()) * 1e6,
                                       amplitude_vpp=float(f["amplitude"].get()), phase_deg=float(f["phase"].get()))
            settings.validate(self.config["channels"][name])
        except (ValueError, KeyError) as e:
            messagebox.showerror("Invalid Live Settings (Nothing Sent)", str(e))
            return "break"
        self.submit(name + " Apply Live", lambda: {name: self.require_rig().apply_live(name, settings)}, self.show_states)
        return "break"

    def output(self, name, enabled):
        def action():
            rig = self.require_rig()
            return {name: rig.set_output(name, enabled, return_state=True)}
        self.submit(name + (" Output ON" if enabled else " Output OFF"), action, self.show_states)

    def stop(self):
        if self.rig is None:
            self.note("Disconnected: unable to send output OFF commands.")
            return
        self.rig.request_stop()
        # Queued behind at most the current operation; stop latch prevents ON.
        def action():
            failures = self.rig.all_off()
            if failures:
                raise RuntimeError("Could not verify these outputs are OFF. Check the instruments:\n" + "\n".join(failures))
        future = self.executor.submit(action)
        future.add_done_callback(lambda f: self.events.put(("All Outputs OFF", f, lambda _: [s.set("Output: OFF; apply settings before enabling") for s in self.states.values()])))
        self.note("Stop requested; waiting for USB confirmation.")

    def align(self, alias):
        self.submit(alias + " Align CH1/CH2", lambda: self.require_rig().align_phase(alias))

    def trigger(self, alias):
        self.submit(alias + " Device Trigger", lambda: self.require_rig().trigger(alias))

    def save_profile(self):
        try:
            value = {n: asdict(self.collect(n)) for n in self.fields}
        except Exception as e:
            messagebox.showerror("Save Failed", str(e))
            return
        filename = filedialog.asksaveasfilename(title="Save Settings", initialdir=BASE, defaultextension=".json", filetypes=[("Settings JSON", "*.json")])
        if filename:
            Path(filename).write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")
            self.note("Settings saved. Output ON states are not saved.")

    def load_profile(self):
        filename = filedialog.askopenfilename(title="Load Settings", initialdir=BASE, filetypes=[("Settings JSON", "*.json")])
        if not filename:
            return
        try:
            data = json.loads(Path(filename).read_text(encoding="utf-8-sig"))
            validated = {}
            for name, values in data.items():
                if name not in self.fields:
                    raise ValueError("Unknown channel: " + name)
                setting = ChannelSettings(**values)
                setting.validate(self.config["channels"][name])
                validated[name] = setting
            for name, s in validated.items():
                values = dict(frequency=s.frequency_hz/1e6, amplitude=s.amplitude_vpp, phase=s.phase_deg,
                              mode=next(k for k, v in MODES.items() if v == s.mode), start=s.sweep_start_hz/1e6,
                              stop=s.sweep_stop_hz/1e6, time=s.sweep_time_s*1000, spacing=s.sweep_spacing,
                              cycles=s.burst_cycles, source=next(k for k, v in SOURCES.items() if v == s.trigger_source))
                for key, value in values.items():
                    self.fields[name][key].set(str(value))
            self.note("Settings loaded into the interface. Nothing has been sent to the instruments.")
        except Exception as e:
            messagebox.showerror("Load Failed", str(e))

    def close(self):
        if self.closing:
            return
        self.closing = True
        if self.rig:
            self.rig.request_stop()
        if not self.busy:
            self.finish_close()
        else:
            self.note("Waiting for the current operation, then turning outputs OFF and exiting.")

    def finish_close(self):
        def action():
            if self.rig:
                self.rig.close()
        future = self.executor.submit(action)
        def wait():
            if not future.done():
                self.root.after(100, wait)
                return
            try:
                future.result()
            except Exception as e:
                messagebox.showerror("Output Shutdown Not Verified", str(e) + "\nCheck the generators or switch off the amplifiers.")
            self.executor.shutdown(wait=False)
            self.root.destroy()
        wait()


def main():
    logging.basicConfig(filename=BASE / "control.log", encoding="utf-8", level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    configure_dpi()
    root = tk.Tk()
    try:
        App(root)
    except Exception as e:
        messagebox.showerror("Startup Failed", str(e))
        root.destroy()
        return
    root.mainloop()


if __name__ == "__main__":
    main()
