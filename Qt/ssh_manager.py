import paramiko
from io import StringIO
import json

class SSHManager:
    """网页8/9/10的Paramiko实现"""
    def __init__(self, host, user, pwd):
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.ssh.connect(host, username=user, password=pwd)
        self.sftp = self.ssh.open_sftp()
    
    def get_json(self, remote_path):
        """下载JSON文件"""
        with self.sftp.file(remote_path, 'r') as f:
            return json.load(f)
    
    def put_json(self, data, remote_path):
        """上传JSON文件"""
        buffer = StringIO()
        json.dump(data, buffer, indent=4)
        self.sftp.putfo(buffer, remote_path)
    
    def exec_command(self, cmd):
        """执行Linux命令"""
        stdin, stdout, stderr = self.ssh.exec_command(cmd)
        return stdout.read().decode()