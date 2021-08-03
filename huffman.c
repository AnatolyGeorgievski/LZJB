#include <stdint.h>
/*!
https://www.compression.ru/download/articles/huff/huffman_1952_minimum-redundancy-codes.pdf
https://webspace.science.uu.nl/~leeuw112/huffman.pdf
http://compression.ru/compression.ru/download/huff.html

Порядок вычисления кодов.

Рассчитать распределение кодов, частоту использования.
Сформировать вектор весовых коэффициентов по массиву упорядоченному по убыванию
*/
struct _Smap {
    uint32_t dist:12;// смещение
    uint32_t mlen:11;// длина кода
    uint32_t nlit:9;// длина литерала
};
/*!
    \param a - алфавит кодовых слов вывернутый по частоте использования,
    \param code_freq -- частота использования по кодам
    \param code_length - длина бит на код
    \param cl_len - длина алфавита, число используемых кодов.
    bl_count -- распределение кодов по длинам бит (количество кодов в алфавите, которые имеют длину i
    bl_max -- номер кода.
 */


typedef struct _HTreeNode HTreeNode_t;
struct _HTreeNode {
    uint32_t weight;
    uint16_t left;
    uint16_t right;
};
typedef struct _HQueue HQueue_t;
struct _HQueue {
    uint32_t weight;
    uint16_t node_id;
    uint16_t prev;
};
#define MAX_BITS 15
#define UNDEF 0xFFFF
/*! \brief распределение по длинам */
static void code_length_count(uint8_t *code_lengths, int cl_len, uint8_t *bl_count)
{
    int i;
    for (i = 0; i <= MAX_BITS; i++)
        bl_count[i]=0;
    for (i = 0; i < cl_len; i++) {
        bl_count[code_lengths[i]]++;
    }
}
/*! \brief оценка степени сжатия

    Итоговую величину следует поделить на число бит исходных кодов (8)

 */
float huffman_estimate(const uint8_t *code_lengths, const uint16_t *weights, int cl_len)
{
    uint32_t sum=0,sum1=0;
    int i;
    for(i=0;i<cl_len; i++) {
        sum1+= weights[i];
        sum += weights[i]*code_lengths[i];
    }
    return (float)sum/sum1;
}
/*! \brief составляет таблицу кодирования из таблицы длины кодов и длин */
void huffman_gen_codes(uint8_t *code_lengths, int cl_len, uint8_t *bl_count, uint16_t *codes)
{// построение таблицы алфавита из max_code
    uint16_t code = 0;
    uint16_t next_code[MAX_BITS+1];
    next_code[0]=0;
    bl_count [0]=0;
    int i;
    for (i = 1; i < MAX_BITS; i++) {
        next_code[i] = code;
        code = (code + bl_count[i]) << 1;
    }
    next_code[i] = code;

    for (i=0; i< cl_len; i++) {
        int len = code_lengths[i];
        if (len != 0) {
            codes[i] = next_code[len]++;
        }
    }
}
/*! \brief строит дерево и заполняет таблицы длины кодов

Сами коды не генерируются, сами коды в отдельном расчете. см. \b huffman_gen_codes
    \param code_lengths - таблица длин кодов, заполняется в алгоритме
    \param weights - весовые коэффициенты - частота использования кода
    \param cl_len длина (размер) кодовой таблицы, в кодовой таблице
    \return длины кодов в массиве code_lengths.
 */
void huffman_tree(uint8_t *code_lengths, const uint16_t *weights, int cl_len)
{
    uint16_t tree_alloc_id=0;// номер для выделения структуры в массиве.
    HTreeNode_t nodes[cl_len+cl_len];// сумма геометрической прогрессии для 2*b(1-1/2^n) округлил до 2*b
    // nested functions
    uint16_t tree_node_new (uint16_t code, uint16_t right_code, uint32_t weight) {
        HTreeNode_t *node = &nodes[tree_alloc_id];
        node->left=code;
        node->right=right_code;// для конечных нод справа будет символ UNDEF
        node->weight=weight;
        return tree_alloc_id++;
    }
    // коды можно генерировать из обхода этого дерева.
    void tree_gen_codes(uint16_t node_id, int code_len){
        uint16_t right = nodes[node_id].right;
        uint16_t left  = nodes[node_id].left;
        if (right==UNDEF){
            code_lengths[left] = code_len;
            return;
        }
        code_len++;
        tree_gen_codes(left,  code_len);
        tree_gen_codes(right, code_len);
    }

    HQueue_t queue[cl_len];
    int queue_len=0;
    uint16_t queue_alloc_id=0;
    uint16_t top = UNDEF;// вершина стека
    uint16_t priority_queue_pop ()
    {
        queue_len--;
        queue_alloc_id = top;// освобождаемый элемент очереди
        top = queue[top].prev;
        return queue[queue_alloc_id].node_id;
    }
    void priority_queue_push(uint16_t node_id, uint32_t weight)
    {
        uint16_t id = queue_alloc_id++;
        queue[id].node_id = node_id;
        queue[id].prev = top;
        queue[id].weight = weight;
        queue_len++;

        uint16_t prev = top, next=UNDEF;
        while (prev!=UNDEF && queue[prev].weight < weight) {
            next = prev;
            prev = queue[prev].prev;
        }
        queue[id].prev = prev;
        if (next != UNDEF)
            queue[next].prev = id;
        else
            top = id;
    }
    int i;
    for(i=0;i< cl_len; i++) {// по всем символам ненулевым
        //uint16_t n = a[i];
        uint32_t weight = weights[i];
        uint16_t node_id = tree_node_new(i, UNDEF, weight);
        priority_queue_push (node_id, weight);
    }
    uint16_t left, right;
    while(1) {
        left  = priority_queue_pop();
        right = priority_queue_pop();
        if (queue_len==0) break;
        uint32_t weight = nodes[left].weight + nodes[right].weight;
        uint16_t node_id  = tree_node_new(left, right, weight);
        priority_queue_push (node_id, weight);
    }
    tree_gen_codes(left,  1);
    tree_gen_codes(right, 1);
}


/*! \brief выполнить кодирование сегмента */
uint8_t* huffman_fixed_encode(uint8_t *dst, uint8_t *src, struct _Smap *map, int map_len)
{
    int i, n_bits=3;
    int bits, extra, index;
    uint32_t code, nlit;
    uint32_t stream=3;// btype|final=01|1
    // nested functions
    void stream_add_bits(uint32_t code, int bits){
        stream |= (code&((1<<bits)-1))<<n_bits;
        n_bits += bits;
        while (n_bits>=8) {
            *dst++ = (uint8_t)stream;
            stream>>=8;
            n_bits -=8;
        }
    }
    void huffman_set_bits(uint32_t code, int bits){
        do {
            stream |= (code&1) << n_bits++;
            code>>=1;
        } while(--bits);
        while (n_bits>=8) {
            *dst++ = (uint8_t)stream;
            stream>>=8;
            n_bits -=8;
        }
    }
    for(i=0; i<map_len; i++){
        nlit = map[i].nlit;
        if (nlit) do{
            code = *src++;
            if (code<144) {
                code = code + 0b00110000;
                huffman_set_bits(code, 8);
            } else {
                code = code-144 + 0b110010000;
                huffman_set_bits(code, 9);
            }
        } while (--nlit);
        code = map[i].mlen;
        if (code==0) {
            huffman_set_bits(0, 7);
            break;
        }
        src += code;// пропуск
        code -= 2;//3; у нас 2 байта, в Deflate 3 байта.
        if (0 && code>=255) { // если длина >255 то extra=0, index=28, выполнить разбиение!!!
            extra = 0;
            index = 28;
            code  = 255;
            // todo вставить разбиение
            //continue;
        }
        bits  = code==0?0:32-__builtin_clz(code);
        extra = bits>2?bits-2-1:0;
        index = (extra<<2)+(code>>extra); // индекс по таблице
        if (index<(280-256))
            huffman_set_bits(index, 7);
        else
            huffman_set_bits(index-(280-256)+0b11000000, 8);
        if (extra) {
            stream_add_bits(code, extra);
        }
        code = map[i].dist;
        bits  = code==0?0: 32-__builtin_clz(code);
        extra = bits>1?bits-1-1:0;
        index = (extra<<1)+(code>>extra);// код по таблице
        huffman_set_bits(index, 5);
        if (extra) {
            stream_add_bits(code, extra);
        }
        // stream_flush
    }
    if (n_bits){// выгрузить остаток потока
        while(n_bits>0) {
            *dst++ = stream;
            stream>>=8;
            n_bits -=8;
        }
    }
    return dst;
}
typedef struct _huffman huffman_t;
struct _huffman {
    uint16_t* codes;
    uint8_t * code_lengths;
    uint16_t* dl_codes;
    uint8_t * dl_code_lengths;
};
/*! \brief выполнить кодирование сегмента */
uint8_t* deflate_encode(uint8_t *dst, uint8_t *src, struct _Smap *map, int map_len, huffman_t* ctx)
{
    int i, n_bits=3;
    int bits, extra, index;
    uint32_t code, nlit;
    uint32_t stream=4;// btype|final=10|1
    // nested functions
    void stream_add_bits(uint32_t code, int bits){
        stream |= (code&((1<<bits)-1))<<n_bits;
        n_bits += bits;
        while (n_bits>=8) {
            *dst++ = (uint8_t)stream;
            stream>>=8;
            n_bits -=8;
        }
    }
    void huffman_encode(int n, uint16_t* codes, uint8_t* code_lens){
        uint32_t code = codes[n];// коды предварительно развернуть
        int bits = code_lens[n];
        stream_add_bits(code, bits);
    }

    for(i=0; i<map_len; i++){
        nlit = map[i].nlit;
        if (nlit) do{
            index= *src++;
            huffman_encode(index, ctx->codes, ctx->code_lengths);
        } while (--nlit);
        code = map[i].mlen;
        if (code==0) {
            huffman_encode(256, ctx->codes, ctx->code_lengths);
            break;
        }
        src += code;// пропуск
        code -= 3;// у нас 2 байта, в Deflate 3 байта.
        bits  = code==0?0:32-__builtin_clz(code);
        extra = bits>=3?bits-3:0;
        index = (extra<<2)+(code>>extra); // индекс по таблице
        huffman_encode(index+257, ctx->codes, ctx->code_lengths);
        if (extra) {
            stream_add_bits(code, extra);
        }
        code = map[i].dist;
        bits  = code==0?0: 32-__builtin_clz(code);
        extra = bits>=2?bits-2:0;
        index = (extra<<1)+(code>>extra);// код по таблице
        huffman_encode(index, ctx->dl_codes, ctx->dl_code_lengths);
        if (extra) {
            stream_add_bits(code, extra);
        }
        // stream_flush
    }
    if (n_bits){// выгрузить остаток потока
        while(n_bits>0) {
            *dst++ = stream;
            stream>>=8;
            n_bits -=8;
        }
    }
    return dst;
}
#if defined(TEST_HUFFMAN)
#include <stdio.h>
int main()
{
    uint16_t weights[13]={20,18,10,10, 10,6,6,4, 4,4,4,3,1};
    uint16_t codes[16];
    uint8_t  code_lengths[16];
    int cl_len = 13;
    huffman_tree(code_lengths, weights, cl_len);
    float ratio =
    huffman_estimate(code_lengths, weights, cl_len);

    uint8_t bl_count[MAX_BITS+1];
    code_length_count(code_lengths, cl_len, bl_count);
    huffman_gen_codes(code_lengths, cl_len, bl_count, codes);
    printf("No |wght|bits|code\n");
    int i,j;
    char code[MAX_BITS+1];
    for(i=0;i<cl_len; i++) {
        for(j=0;j<code_lengths[i]; j++)
            code[j] = codes[i]&(1<<(code_lengths[i]-j-1))?'1':'0';
        printf("%2d |%3d |%3d |%-.*s\n", i, weights[i], code_lengths[i], code_lengths[i], code);
    }
    printf("ratio = %1.2f%%\n", ratio/4);

    uint16_t weights2[10]={20,17,6,3, 2,2,2,1, 1,1,};
    cl_len = 10;
    huffman_tree(code_lengths, weights2, cl_len);
    ratio =
    huffman_estimate(code_lengths, weights2, cl_len);

    code_length_count(code_lengths, cl_len, bl_count);
    huffman_gen_codes(code_lengths, cl_len, bl_count, codes);
    printf("No |wght|bits|code\n");
    for(i=0;i<cl_len; i++) {
        for(j=0;j<code_lengths[i]; j++)
            code[j] = codes[i]&(1<<(code_lengths[i]-j-1))?'1':'0';
        printf("%2d |%3d |%3d |%-.*s\n", i, weights2[i], code_lengths[i], code_lengths[i], code);
    }
    printf("ratio = %1.2f%%\n", ratio/4);

    //uint16_t weights3[7]={369,245,100,42, 24,19,6};
    uint16_t weights3[7]={369,245,100,12, 9,7,6};
    cl_len = 7;
    huffman_tree(code_lengths, weights3, cl_len);
    ratio =
    huffman_estimate(code_lengths, weights3, cl_len);

    code_length_count(code_lengths, cl_len, bl_count);
    huffman_gen_codes(code_lengths, cl_len, bl_count, codes);
    printf("No |wght|bits|code\n");
    for(i=0;i<cl_len; i++) {
        for(j=0;j<code_lengths[i]; j++)
            code[j] = codes[i]&(1<<(code_lengths[i]-j-1))?'1':'0';
        printf("%2d |%3d |%3d |%-.*s\n", i, weights3[i], code_lengths[i], code_lengths[i], code);
    }
    printf("ratio = %1.2f%%\n", ratio/3);

    return 0;
}

#endif // defined
