#include <stdio.h>

int main(int argc, char** argv) {
    FILE* f = fopen(argv[1], "r");
    char c = 0, last = 1;
    int cnt = 0;
    while(c != EOF) {
        c = fgetc(f);
        if (c != last) {
            cnt++;
            last = c;
        }
    }
    printf("%d\n", cnt);
}
