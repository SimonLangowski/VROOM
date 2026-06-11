class PrecomputeIO:
    def __init__(self, filename):
        self.filename = filename
        self.file = open(filename, "w")
        self.write_counter = 0

    def vector(self, vec):
        for v in vec:
            print(v, file=self.file, end=" ")
            self.write_counter += 1
        print("", file=self.file)

    def matrix(self, mat):
        for row in mat:
            self.vector(row)
    
    def scalar(self, s):
        self.vector([s])

