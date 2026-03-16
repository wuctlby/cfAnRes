import subprocess

class PdgEnv:
    def __init__(self, batch_mode=False):
        """Initialize, start bash and load environment"""
        self.proc = subprocess.Popen(
            "/bin/bash", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        self.stdin = self.proc.stdin
        self.stdout = self.proc.stdout
        self.batch_mode = batch_mode

        # Load environment
        self._load_environment()

        # Discard environment loading output until 'Environment loaded'
        self._discard_env_output()

    def _load_environment(self):
        """Send commands to load environment"""
        self.stdin.write("source ~/.bashrc\n")
        self.stdin.write("alienv enter O2/latest O2DPG/latest\n")
        self.stdin.write("echo 'Environment loaded'\n")
        self.stdin.flush()

    def _discard_env_output(self):
        """Discard environment loading output until 'Environment loaded' appears"""
        while True:
            line = self.stdout.readline()
            if not line:
                break
            if "Environment loaded" in line:
                break
        if not self.batch_mode:
            print("Environment loaded, ready to run commands.")

    def _print_output(self, output):
        """Print command output"""
        output_lines = []
        while True:
            line = output.readline()
            if not line:
                break
            if line.strip() == "---COMMAND END---":
                break
            if 'O2DPG' in line:
                continue
            output_lines.append(line)
            if not self.batch_mode:
                print(line, end="")
        return "".join(output_lines)

    def run(self, command):
        """Run a command and return its output"""
        self.stdin.write(command + "\n")
        self.stdin.write("echo '---COMMAND END---'\n")
        self.stdin.flush()
        return self._print_output(self.stdout)

    def close(self):
        """Close bash process"""
        if self.stdin:
            self.stdin.write("exit\n")
            self.stdin.flush()
            self.stdin.close()
        if self.proc:
            self.proc.wait()