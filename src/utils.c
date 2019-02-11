void mod3(int input, int res[2]){

    switch (input)
    {
        case 0:
            res[0] = 1;
            res[1] = 2;
            break;
        case 1:
            res[0] = 2;
            res[1] = 0;
            break;
        case 2:
            res[0] = 0;
            res[1] = 1;
            break;
        default:
            break;
    }

}

double norm1d(double r[], int n){
    double sum = 0.0;
    int i;
    for (i = 0; i < n; i++){
        sum = sum + r[i] * r[i];
    }
    return sqrt(sum);
}

double norm2d(double r[][4], int n){
    double sum = 0.0;
    int i, j;
    double row[4];
    for (i = 0; i < n; i++){
        for (j = 0; j < 4; j++){
            row[j] = r[i][j];
        } 
        sum = sum + norm1d(row, 4) * norm1d(row, 4);
    }
    return sqrt(sum);
}

