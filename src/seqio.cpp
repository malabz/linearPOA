#include "seqio.hpp"

void gap_insert(char *raw_str, char *ans_str, size_t raw_len, size_t ans_len)
{
    size_t required_raw_len = 0;
    for (size_t i = 0; i < ans_len; i++) {
        if (ans_str[i] != '-') required_raw_len++;
    }
    if (required_raw_len != raw_len) {
        fprintf(stderr, "Warning: ans_str requires %zu non-dash chars, but raw_len=%zu. May generate uncorrect answer.\n", required_raw_len, raw_len);
    }

    size_t raw_i = 0, ans_i = 0;
    for(; raw_i < raw_len && ans_i < ans_len; )
    {
        if(ans_str[ans_i] == '-') { ++ ans_i; continue; }
        ans_str[ans_i] = raw_str[raw_i];
        ++ ans_i; ++ raw_i;
    }
    // fprintf(stderr, "raw_i = %ld, raw_len = %ld, ans_i = %ld, ans_len = %ld\n", raw_i, raw_len, ans_i, ans_len);
    // assert(! raw_len);
}
