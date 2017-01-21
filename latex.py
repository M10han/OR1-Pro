import numpy
def start_latex():
#Beginning the latex initialization
    latex = (r"\documentclass{article}"
             r"\usepackage{amsmath}"
             r"\begin{document}"
             r"\title{Revised Simplex Method Solution}"
             r"\maketitle"
             r"\pagestyle{plain}")
    return latex


def print_objective_function(A,b,c,latex):

    A = numpy.array(A)
    b = numpy.array(b)
    c = numpy.array(c)
    latex += (r"\section{Standard Form}")
    latex += r'Maximize: $ '
    temp_text=''
    for i in range(len(c)):
        temp_text += (str(int(c[i])) + "x_" + str(i + 1) + "+")
    latex += temp_text[0:len(temp_text) - 1]
    latex+=r"$"

    latex += r"\newline Objective Function Constraints:\\ \newline $"

    temp_text=""
    for i in range(len(b)):
        b_value = b[i]
        for j in range(len(A[i])):
            temp_text += (str(int(A[i].item(j)))+r"x_"+str(j+1)+"+")
        latex += temp_text[0:len(temp_text)-1]+r"<="+str(int(b_value))+r"\newline"
        temp_text=""
    latex+=r"$ \newline"
    latex+=r"Printing matrices :\newline \\"
    latex += r"\["
    latex += r"A="
    latex += r"\begin{bmatrix} "
    A_row = ""
    for ele in A:
        for i in range(len(ele)):
            A_row += (str(ele[i]) + r" & ")
        A_row = A_row[0:len(A_row) - 2]
        A_row = A_row + r" \\"
        latex += A_row
        A_row = ""
    latex += r"\end{bmatrix}"
    latex += r"\]"
    latex += r"\["
    latex += r"b="
    latex += r"\begin{bmatrix} "
    b_row = ""
    for ele in b:
        for i in range(len(ele)):
            b_row += (str(ele[i]) + r" & ")
        b_row = b_row[0:len(b_row) - 2]
        latex += b_row
        b_row = ""
    latex += r"\end{bmatrix}"
    latex += r"\]"
    latex += r"\["
    latex += r"c="
    latex += r"\begin{bmatrix} "
    C_row = ""
    for ele in c:
        for i in range(len(ele)):
            C_row += (str(ele[i]) + r" & ")
        C_row = C_row[0:len(C_row) - 2]
        C_row = C_row + r" \\"
        latex += C_row
        C_row = ""
    latex += r"\end{bmatrix}"
    latex += r"\]"

    return latex

def print_bnstar_matrix(mat,latex,tex):

    mat=numpy.array(mat)
    latex += r"$ \newline"

    latex += r"\["
    latex += (tex+r"=")
    latex += r"\begin{bmatrix} "
    row = ""
    for ele in mat:
        for i in range(len(ele)):
            row += (str(ele[i]) + r" & ")
        row = row[0:len(row) - 2]
        row = row + r" \\"
        latex += row
        row = ""
    latex += r"\end{bmatrix}"
    latex += r"\]"

    return latex

def save_latex(file_path,latex):
    f=open(file_path,"w+")
    f.write(latex)
    f.close()