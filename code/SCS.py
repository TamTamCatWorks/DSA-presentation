import heapq

MOD = 10**9 + 7
BASE = 311


class StringHash:
    """Cấu trúc Rolling Hash để lấy mã băm đoạn con trong O(1)"""

    def __init__(self, s=""):
        n = len(s)
        self.hash = [0] * (n + 1)
        self.power = [1] * (n + 1)

        for i in range(n):
            self.power[i + 1] = (self.power[i] * BASE) % MOD
            self.hash[i + 1] = (self.hash[i] * BASE + ord(s[i])) % MOD

    def get_hash(self, left, right):
        """Lấy hash của đoạn [left, right] (tính từ 1)"""
        length = right - left + 1
        val = (self.hash[right] - self.hash[left - 1] * self.power[length]) % MOD
        return (val + MOD) % MOD


class MergeResult:
    """Kết quả khi ghép 2 chuỗi"""

    def __init__(self, overlap=0, merged=""):
        self.overlap = overlap
        self.merged = merged


def find_best_merge(s1, s2):
    """Tìm phần trùng dài nhất giữa hậu tố của s1 và tiền tố của s2"""
    n1 = len(s1)
    n2 = len(s2)
    max_len = min(n1, n2)

    h1 = StringHash(s1)
    h2 = StringHash(s2)

    for length in range(max_len, 0, -1):
        suffix_hash = h1.get_hash(n1 - length + 1, n1)
        prefix_hash = h2.get_hash(1, length)

        if suffix_hash == prefix_hash and s1[n1 - length:] == s2[:length]:
            return MergeResult(length, s1 + s2[length:])

    return MergeResult(0, s1 + s2)



def remove_substrings(strings):
    """Loại bỏ các chuỗi nằm gọn hoàn toàn bên trong chuỗi khác"""
    n = len(strings)
    removed = [False] * n

    for i in range(n):
        if removed[i]:
            continue

        for j in range(n):
            if i == j or removed[j]:
                continue

            if strings[i] in strings[j]:
                # Nếu 2 chuỗi giống hệt nhau thì chỉ xóa chuỗi đứng sau
                if strings[i] == strings[j]:
                    removed[max(i, j)] = True
                else:
                    removed[i] = True
                    break

    return [strings[i] for i in range(n) if not removed[i]]



def push_heap(heap, matrix, versions, i, j):
    """Đẩy vào heap"""
    overlap = matrix[i][j].overlap
    heapq.heappush(heap, (-overlap, i, j, versions[i], versions[j]))


def init_merge_data(strings, versions):
    """Tạo ma trận kết quả ghép ban đầu và đưa tất cả các cặp vào heap"""
    n = len(strings)
    matrix = [[MergeResult() for _ in range(n)] for _ in range(n)]
    heap = []

    for i in range(n):
        for j in range(n):
            if i != j:
                matrix[i][j] = find_best_merge(strings[i], strings[j])
                push_heap(heap, matrix, versions, i, j)

    return matrix, heap


def pop_best_pair(heap, used, versions):
    """Lấy cặp (i, j) có độ trùng lớn nhất và còn hợp lệ ra khỏi heap"""
    while True:
        _, i, j, ver_i, ver_j = heapq.heappop(heap)

        # Bỏ qua nếu 1 trong 2 chuỗi đã bị gộp mất
        if used[i] or used[j]:
            continue
        # Bỏ qua nếu phiên bản lưu trong heap đã cũ
        if ver_i != versions[i] or ver_j != versions[j]:
            continue

        return i, j


def update_after_merge(u, strings, matrix, heap, used, versions):
    """Tính lại độ trùng của chuỗi u (vừa mới gộp) với các chuỗi v còn lại"""
    n = len(strings)
    for v in range(n):
        if used[v] or v == u:
            continue

        # Cập nhật cả 2 chiều: u -> v và v -> u
        matrix[u][v] = find_best_merge(strings[u], strings[v])
        matrix[v][u] = find_best_merge(strings[v], strings[u])

        push_heap(heap, matrix, versions, u, v)
        push_heap(heap, matrix, versions, v, u)


def build_shortest_superstring(strings):
    """Thuật toán chính: Ghép dần các cặp có độ trùng lớn nhất"""
    if not strings:
        return ""

    strings = remove_substrings(strings)
    n = len(strings)
    if n == 0:
        return ""
    if n == 1:
        return strings[0]

    # used: đánh dấu chuỗi đã bị gộp vào chuỗi khác
    # versions: lưu phiên bản của chuỗi để lọc bỏ các cặp cũ trong heap
    ## Vì trong list string sau khi gộp, tại vị trí i sẽ có chuỗi mới
    ### Nên các cặp ghép cũ với chuỗi i sẽ không còn hợp lệ nữa.
    #### Mà trong heap không thể đánh dấu trực tiếp các cặp cũ, nên ta dùng versions để kiểm tra.
    used = [False] * n
    versions = [0] * n
    matrix, heap = init_merge_data(strings, versions)

    # Có n chuỗi thì cần đúng n - 1 bước gộp
    for _ in range(n - 1):
        i, j = pop_best_pair(heap, used, versions)

        # Gộp chuỗi j vào chuỗi i
        strings[i] = matrix[i][j].merged
        used[j] = True
        versions[i] += 1

        # Cập nhật lại các cặp ghép với chuỗi i mới
        update_after_merge(i, strings, matrix, heap, used, versions)

    # Trả về chuỗi duy nhất chưa bị gộp
    for i in range(n):
        if not used[i]:
            return strings[i]

    return ""


def input_count():
    """Nhập số lượng chuỗi (1 -> 100)"""
    print("Nhap so luong chuoi, khong qua 100: ", end="")
    while True:
        try:
            n = int(input())
            if 1 <= n <= 100:
                return n
        except ValueError:
            pass
        print("So luong chuoi phai nam trong khoang 1 den 100. Nhap lai: ", end="")


def input_single_string():
    """Nhập và kiểm tra 1 chuỗi (độ dài 1 -> 50)"""
    while True:
        s = input().strip()
        if 1 <= len(s) <= 50:
            return s
        print("Chuoi phai co do dai tu 1 den 50. Nhap lai: ", end="")


def input_string_list(n):
    """Nhập danh sách n chuỗi"""
    print("Nhap cac chuoi, moi chuoi dai khong qua 50 ky tu:")
    strings = []
    for _ in range(n):
        strings.append(input_single_string())
    return strings


def main():
    n = input_count()
    strings = input_string_list(n)

    result = build_shortest_superstring(strings)

    print(f"\nShortest Superstring: {result}")
    print(f"Do dai: {len(result)}")


if __name__ == "__main__":
    main()