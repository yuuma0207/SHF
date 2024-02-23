import subprocess
import re
import os
import time
import shlex

def run_plotall_script(folder_name):
    command = ["python3", "./../src/plotall.py", folder_name]
    subprocess.Popen(command)
    return None


def extract_folder_name(text):
    # フォルダ名のパターン（例: 0124/113539/）を探します
    match = re.search(r'\d{4}/\d{6}/', text)
    if match:
        return match.group(0)  # マッチしたフォルダ名を返します
    return None


def main():
    process = subprocess.Popen(['./main'], stdout=subprocess.PIPE, text=True)
    folder_name = None
    file = None
    not_word = ["plot_wf", "plot_density", "plot_d_wf", "plot_potential", "plot_energy"]
    plot_word = "python3 ./../src/plotall.py"

    while True:
        output = process.stdout.readline()
        # 終了条件
        rc = "endrecord"
        if(output.strip() == rc):
            break

        if output:
            if plot_word in output:
                args = shlex.split(output)
                script_name = args[1]
                folder_name = args[2]
                script_args = args[3:]
                if folder_name:
                    folder_name = os.path.join("./../res", folder_name)
                    run_command = [args[0], script_name, folder_name] + script_args
                    subprocess.Popen(run_command)
            # フォルダ名を含むかチェック
            if folder_name is None:
                folder_name = extract_folder_name(output.strip())

                # フォルダ名が見つかったらファイルを開く
                if folder_name:
                    base_path = "./../res"
                    folder_name = os.path.join(base_path, folder_name)
                    file_path = os.path.join(folder_name, 'output.txt')
                    file = open(file_path, 'w')
                    if file:
                        with open('./../src/input_parameters.txt', 'r') as input_file:  # 1. 読み込みたいファイルを開く
                            content_to_add = input_file.read()  # 2. 内容を読み込む
                        file.write(content_to_add)  # 4. 読み込んだ内容を追記する
                    else:
                        print("ファイルが開けませんでした")
                        exit(1)

            # ファイルが開かれていれば書き込む
            if file:
                if not any(word in output for word in not_word):
                    file.write(output)
                    file.flush()

    # ファイルがまだ開かれている場合は閉じる
    if file and not file.closed:
        file.close()


if __name__ == "__main__":
    main()