import datetime


class Logger:
    def __init__(self, PATH, filename, erase_on_init=False):
        if PATH == "":
            self.output_file_path = filename
        else:
            self.output_file_path = f"{PATH}/{filename}"

        # Open the file in write mode, which will erase it
        if erase_on_init:
            with open(self.output_file_path, "w") as file:
                pass # Do nothing so the file is empty
        return

    def log(self, msg):
        """Write a message to the log file"""
        current_time = datetime.datetime.now()
        with open(self.output_file_path, "a") as file:
            file.write(current_time.strftime("%Y/%m/%d %H:%M:%S : ") + msg + "\n")
        return