# convert .msh grid file into GridNodes.dat and GridElements.dat pair
def format_mesh(mesh_file, nodes_file, elements_file):
    with open(mesh_file, "r") as mf, open(nodes_file, "w") as nf, open(elements_file, "w") as ef:
        regime = 0
        number_of_elements = 0
        for line in mf:
            if line.startswith("$Nodes"):
                regime = 1
            elif line.startswith("$Elements"):
                regime = 2
            
            if line[0] == "$":
                continue
            
            if regime == 1:
                nf.write(line)
            elif regime == 2:
                lineElements = line.split()
                if len(lineElements) == 9:
                    ef.write(" ".join(lineElements[-4:]) + "\n")
                    number_of_elements += 1
                    
    with open(elements_file, 'r+') as ef:
	    lines = ef.readlines()
	    ef.seek(0)
	    ef.writelines( [str(number_of_elements)+"\n"]+lines )
                
if __name__ == "__main__":
    format_mesh("disk.msh", "GridNodes.dat", "GridElements.dat")

