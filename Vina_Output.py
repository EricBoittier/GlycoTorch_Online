class PDBQTOutput:
    def __init__(self, name, model, path):
        self.name = name
        self.model = model.replace(" ", "_")
        self.model = self.model.strip("\n")+"_"
        self.path = path
        self.lines = []


    def write_model(self):

        filename = self.name.split("/")[-1]
        filepath = self.model+filename.strip(" ")

        output = open(self.path+"/"+filepath, "w")

        for line in self.lines:
            output.write(line)

    def add_line(self, line):
        self.lines.append(line)


class OutputFile:
    """
    A class to handle files
    """
    def __init__(self, filename, path="."):
        self.pdbqt_files = {}
        self.filename = filename
        self.file_lines = []
        self.read_file()
        self.path = path

    def set_filename(self, string):
        self.filename = string

    def get_filename(self):
        return self.filename

    def get_file(self):
        return self.file_lines

    def read_file(self):
        with open(self.filename, 'r') as f:
            self.file_lines = f.readlines()

    def make_pdbqt_files(self, location="."):
        filename = self.filename.split("/")[-1]

        pdbqt_outputs = []

        add_model = True
        for line in self.file_lines:

            if line.__contains__("ENDMDL"):
                add_model = True

            if line.__contains__("MODEL"):
                pdbqt_outputs.append(PDBQTOutput(filename, line, location))
                add_model = False

            if not add_model:
                pdbqt_outputs[-1].add_line(line)

        for output in pdbqt_outputs:
            output.write_model()

    def get_pdbqt_files(self):
        return self.pdbqt_files


