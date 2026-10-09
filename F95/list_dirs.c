/* List immediate subdirectory basenames of root into a flat buffer.
 * names_buf: capacity names_cap entries of name_len chars each (not NUL-padded
 *            beyond the name; caller treats as Fortran character(name_len)).
 * Returns number of subdirs written, or -1 on error.
 * Skips "." and "..". Only includes entries that are directories (or symlink-to-dir).
 */
#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

int ossssim_list_subdirs(const char *root, char *names_buf, int names_cap,
                         int name_len)
{
    DIR *dir;
    struct dirent *ent;
    int n = 0;
    char full[4096];
    struct stat st;
    size_t root_len, namelen;

    if (!root || !names_buf || names_cap <= 0 || name_len <= 0)
        return -1;

    dir = opendir(root);
    if (!dir)
        return -1;

    root_len = strlen(root);
    while ((ent = readdir(dir)) != NULL) {
        if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0)
            continue;
        namelen = strlen(ent->d_name);
        if ((int)namelen >= name_len)
            continue;
        if (root_len + 1 + namelen >= sizeof(full))
            continue;
        memcpy(full, root, root_len);
        full[root_len] = '/';
        memcpy(full + root_len + 1, ent->d_name, namelen + 1);
        if (stat(full, &st) != 0)
            continue;
        if (!S_ISDIR(st.st_mode))
            continue;
        if (n >= names_cap) {
            closedir(dir);
            return -1;
        }
        memset(names_buf + n * name_len, ' ', (size_t)name_len);
        memcpy(names_buf + n * name_len, ent->d_name, namelen);
        n++;
    }
    closedir(dir);
    return n;
}

int ossssim_file_exists(const char *path)
{
    struct stat st;
    if (!path)
        return 0;
    if (stat(path, &st) != 0)
        return 0;
    return S_ISREG(st.st_mode) ? 1 : 0;
}
