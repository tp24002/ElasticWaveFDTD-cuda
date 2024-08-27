#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// ここでの"../header/struct.h"と"../header/init.h"はユーザーが作成したヘッダーファイルであると仮定しています。
// 必要に応じて、ヘッダーファイルのパスを適切に調整してください。
#include "../header/struct.h"

// 乱数を生成する関数
void random(Coord con_size, Coord *clack, int ratio){
    // printf("ratio:%d%\n",ratio);
    if(con_size.x < 3 || con_size.y < 3 || con_size.z < 3) {
        printf("Cannot place defects.\n");
        return;
    }
    int count = 0;
    // コンクリートセル数
    int max_Patern = con_size.x * con_size.y * con_size.z;
    // 内部欠陥パターン数
    int max_ClackPatern = (con_size.x - 2) * (con_size.y - 2) * (con_size.z - 2);
    // 割合による欠陥数
    int clack_count = max_Patern * ratio / 100;
    if(clack_count > max_ClackPatern){
        printf("The number of internal defects is insufficient.\n");
        return;
    }
    
    // 乱数の種を初期化
    srand(time(NULL));

    while (count < clack_count) {
        // 新しい乱数の組み合わせを生成
        int rand1 = rand() % (con_size.x - 2) + 1;
        int rand2 = rand() % (con_size.y - 2) + 1;
        int rand3 = rand() % (con_size.z - 2) + 1;

        // 重複がないかチェック
        int is_unique = 1;
        for (int i = 0; i < count; i++) {
            if (clack[i].x == rand1 && clack[i].y == rand2 && clack[i].z == rand3) {
                is_unique = 0;
                break;
            }
        }

        // 重複がなければ保存
        if (is_unique) {
            clack[count].x = rand1;
            clack[count].y = rand2;
            clack[count].z = rand3;
            count++;
        }
    }

    
}

int main(){
    int ratio = 10;
    Coord con_size;
    con_size.x = 3;
    con_size.y = 3;
    con_size.z = 3;

    int max_Patern = con_size.x * con_size.y * con_size.z;
    int max_ClackPatern = (con_size.x - 2) * (con_size.y - 2) * (con_size.z - 2);
    int clack_count = max_Patern * ratio / 100;
    // if(clack_count > max_ClackPatern){
    //     printf("aaThe number of internal defects is insufficient.");
    //     return 0;
    // }
    // 構造体の配列を動的に確保する
    Coord *clack = (Coord *)malloc(clack_count * sizeof(Coord));
    if (clack == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    // 乱数を生成
    random(con_size, clack, ratio);
    // 結果を表示
    for (int i = 0; i < clack_count; i++) {
        printf("Pattern %3d: %d, %d, %d\n", i + 1, clack[i].x, clack[i].y, clack[i].z);
    }
    // 使用後はメモリを解放する
    free(clack);

    return 0;
}
