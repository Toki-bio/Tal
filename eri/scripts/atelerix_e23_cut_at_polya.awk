# Cut longA / longB consensus from a plate: start where the original e2-3 row starts, end in the terminal poly-A near column 430.
# Base per column = 50 % consensus row; where it has none, the 35 % row's base if at least half of the copies have a base there.
# The end: first run of >= 8 A (gaps ignored) in the 50 % consensus that begins at column >= 380; kept with TAIL A.
/^>/ { k++; h[k] = substr($0, 2); next }
{ s[k] = s[k] $0 }
END {
    for (i = 1; i <= k; i++) { if (h[i] ~ /35pct/) c35 = i; else if (h[i] ~ /50pct/) c50 = i; else if (h[i] == "e2-3") o = i }
    L = length(s[c50])
    for (j = 1; j <= L; j++) { c = substr(s[o], j, 1); if (c != "-" && c == toupper(c)) { a = j; break } }
    n = 0; q = ""
    for (j = a; j <= L; j++) {
        b = substr(s[c50], j, 1)
        if (b == "-") { have = 0; for (i = 1; i <= k; i++) if (i != c35 && i != c50 && i != o && substr(s[i], j, 1) != "-") have++
                        if (have >= 0.5 * (k - 4)) b = substr(s[c35], j, 1) }
        if (b == "-" || b != toupper(b)) continue
        q = q b; col[++n] = j
    }
    for (p = 1; p <= n - 7; p++) if (col[p] >= 380 && substr(q, p, 8) == "AAAAAAAA") break
    r = p; while (r > 1 && substr(q, r - 1, 1) == "A") r--
    seq = substr(q, 1, r - 1) "AAAAAAAAAA"
    print ">" NAME; print seq
    printf "%s: columns %d (e2-3 start) to poly-A at column %d; %d bp with 10 A\n", NAME, a, col[r], length(seq) > "/dev/stderr"
}
