from PyQt5.QtWidgets import QMainWindow, QApplication
from ui_main import Ui_MainWindow  # 网页1生成的界面类
from ssh_manager import SSHManager  # 网页8的Paramiko实现
import sys
import json

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        # 连接按钮事件（网页3信号槽机制）
        self.ui.btn_load.clicked.connect(self.load_config)
        self.ui.btn_save.clicked.connect(self.save_and_run)
    
    def load_config(self):
        """从服务器加载JSON配置"""
        try:
            self.ssh = SSHManager(
                self.ui.txt_host.text(),
                self.ui.txt_user.text(),
                self.ui.txt_pwd.text()
            )
            data = self.ssh.get_json("/path/server_config.json")  # 网页9文件下载
            self.ui.json_editor.setPlainText(json.dumps(data, indent=4))  # 网页7的JSON处理
        except Exception as e:
            self.ui.log_area.append(f"加载失败: {str(e)}")
    
    def save_and_run(self):
        """保存配置并执行命令"""
        try:
            data = json.loads(self.ui.json_editor.toPlainText())
            self.ssh.put_json(data, "/path/server_config.json")  # 网页9文件上传
            output = self.ssh.exec_command("xxx run")  # 网页8命令执行
            self.ui.log_area.append(f"执行结果:\n{output}")
        except Exception as e:
            self.ui.log_area.append(f"错误: {str(e)}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())