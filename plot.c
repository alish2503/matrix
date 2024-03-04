#include <stdio.h>
#include <stdlib.h>

int main(void) {
    	FILE *gp = popen("gnuplot","w"); 
	if (!gp) {
        	printf("Error opening pipe to GNU plot.\n");
        	exit(0);
    	}
    	fprintf(gp, "set xlabel 'X'\n");
    	fprintf(gp, "set ylabel 'Y'\n");
    	fprintf(gp, "set zlabel 'U'\n");
    	fprintf(gp, "set xrange [0:1]\n");
    	fprintf(gp, "set yrange [0:1]\n");
    	fprintf(gp, "set ticslevel 0\n");
    	fprintf(gp, "set dgrid3d 30,30 qnorm 3\n");
    	fprintf(gp, "set terminal png\n");
    	fprintf(gp, "set output 'heatmap.png'\n");
    	fprintf(gp, "plot 'data.txt' w image\n");
    	fprintf(gp, "set output 'graph.png'\n");
    	fprintf(gp, "set view 40, 40\n");
    	fprintf(gp, "splot 'data.txt' w pm3d notitle\n");
    	pclose(gp);
    	return 0;
}
