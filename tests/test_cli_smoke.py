import subprocess

   def test_cli_help_exits_zero():
       result = subprocess.run(["your-cli", "-h"], capture_output=True)
       assert result.returncode == 0
       assert b"-h" in result.stdout or b"--help" in result.stdout