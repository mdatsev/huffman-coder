#include <getopt.h>
#include <iostream>
#include <filesystem>
#include <sstream>
#include <vector>
#include <queue>
#include <cstring>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>

#include "crc.h"

void set_bit_at(unsigned char &byte, int pos, bool bit)
{
  byte = (byte & ~(1U << (7 - pos))) | (bit << (7 - pos));
}

bool get_bit_at(unsigned char byte, int pos)
{
  return (byte >> (7 - pos)) & 1U;
}

// format:
// dddddddd ... dddddxxx  pppppppp
// d - data
// x - padding
// p - index of first bit of padding in the last byte (1-8), 8 means no padding

class IBitStream
{
public:
  IBitStream(std::istream &stream) : stream(stream)
  {
    // read 2 times to make sure both current and next byte are initialized
    if (!read_byte_from_stream())
    {
      // TODO throw error no padding
    }
    if (!read_byte_from_stream())
    {
      stream_is_empty = true;
    }
  }

  IBitStream(std::istream &stream, uint64_t stream_length) : stream(stream), stream_length(stream_length), has_length(true)
  {
    if (stream_length > 0)
    {
      read_byte_from_stream();
    } else {
      istream_ended = true;
    }
  }

  ~IBitStream()
  {
    //
  }

  bool read_bit()
  {
    if (!has_more())
    {
      // TODO check if this is assert
      return 0;
    }
    bool bit = get_bit_at(has_length ? next_byte : current_byte, next_bit_pos);
    next_bit_pos++;
    if (next_bit_pos == 8 && !istream_ended)
    {
      read_byte_from_stream();
      next_bit_pos = 0;
    }

    curr_pos++;

    return bit;
  }

  bool has_more()
  {
    return has_length
      ? (curr_pos < stream_length)
      : !(stream_is_empty || istream_ended && next_bit_pos >= next_byte);
  }

private:
  bool read_byte_from_stream()
  {
    if (istream_ended)
    {
      // TODO throw, maybe assert
      return false;
    }

    unsigned char read_byte = stream.get();

    if (stream.peek() == EOF)
    {
      istream_ended = true;
    }

    if (stream.bad())
    {
      // TODO throw exception, maybe use this https://gehrcke.de/2011/06/reading-files-in-c-using-ifstream-dealing-correctly-with-badbit-failbit-eofbit-and-perror/
      // for now return
      std::cout << "ERROR 2" << std::endl;
      return false;
    }

    current_byte = next_byte;
    next_byte = read_byte;
    return true;
  }

  unsigned char current_byte = 0;
  unsigned char next_byte = 0; // saving next byte in case it's the padding byte (last byte)
  int next_bit_pos = 0;
  bool stream_is_empty = false;
  bool istream_ended = false;
  std::istream &stream;
  bool has_length = false;
  uint64_t stream_length;
  int curr_pos = 0;
};

class OBitStream
{
public:
  OBitStream(std::ostream &stream) : stream(stream)
  {
  }

  ~OBitStream()
  {
    // TODO assert flushed (padding)
  }

  // flushes remaining bits and writes padding info
  uint64_t flush(bool write_padding = true)
  {
    if (!write_padding) {
      if (next_bit_pos != 0) {
        write_byte_to_stream(current_byte);
      }
      return curr_pos;
    } else {
      if (next_bit_pos == 0)
      {
        write_byte_to_stream(8);
      }
      else
      {
        write_byte_to_stream(current_byte);
        write_byte_to_stream(next_bit_pos);
      }
      return curr_pos;
    }
  }

  void write_bit(bool bit)
  {
    set_bit_at(current_byte, next_bit_pos, bit);

    next_bit_pos++;
    if (next_bit_pos == 8)
    {
      write_byte_to_stream(current_byte);
      next_bit_pos = 0;
    }

    curr_pos++;
  }

  void write_bit_chunk(unsigned char *bits, int length)
  {
    // TODO optimize this
    for (int i = 0; i < length; i++)
    {
      write_bit(get_bit_at(bits[i / 8], i % 8));
    }
  }

private:
  void write_byte_to_stream(unsigned char c)
  {
    stream << c;

    if (stream.fail())
    {
      // TODO throw exception, maybe use this https://gehrcke.de/2011/06/reading-files-in-c-using-ifstream-dealing-correctly-with-badbit-failbit-eofbit-and-perror/
    }
  }

  unsigned char current_byte;
  int next_bit_pos = 0;
  std::ostream &stream;
  int curr_pos = 0;
};

struct HuffmanNode
{
  HuffmanNode(unsigned char byte, int frequency, HuffmanNode *left, HuffmanNode *right)
      : byte{byte}, frequency{frequency}, children{left, right}
  {
  }

  HuffmanNode() {}

  ~HuffmanNode()
  {
    delete children[0];
    delete children[1];
  }
  unsigned char byte;
  int frequency;
  HuffmanNode *children[2] = {nullptr, nullptr};
};

// source https://stackoverflow.com/a/51730733/9156061
void printBT(const std::string &prefix, const HuffmanNode *node, bool isLeft, std::string code)
{
  if (node != nullptr)
  {
    std::cout << prefix;

    std::cout << (isLeft ? "├──" : "└──");

    // print the value of the node
    if (node->children[0] == nullptr && node->children[0] == nullptr)
    {
      std::cout << "[" << node->byte << ":" << node->frequency << "] " << code << std::endl;
    }
    else
    {
      std::cout << "(" << node->frequency << ")" << std::endl;
    }

    // enter the next tree level - left and right branch
    printBT(prefix + (isLeft ? "│   " : "    "), node->children[0], true, code + "0");
    printBT(prefix + (isLeft ? "│   " : "    "), node->children[1], false, code + "1");
  }
}

void printBT(const HuffmanNode *node)
{
  printBT("", node, false, "");
}

namespace Option { 
  unsigned char
    is_raw                 = 0b00000001,
    has_crc                = 0b00000010,
    has_crc_each           = 0b00000100,
    has_compression_levels = 0b00001000,
    has_only_empty_files   = 0b00010000,
    has_frequencies        = 0b00100000,
    version                = 0b00100000;
};

unsigned char global_options = 0;

struct HuffmanBitCode
{
  unsigned char byte;
  int length = 0;
  int capacity = 0;
  unsigned char *bits = nullptr;

  HuffmanBitCode(const HuffmanBitCode &h)
  {
    *this = h;
  }

  HuffmanBitCode &operator=(const HuffmanBitCode &other)
  {
    if (this == &other)
      return *this;

    length = other.length;
    byte = other.byte;
    if (capacity != other.capacity)
    {
      delete[] bits;
      bits = other.capacity ? new unsigned char[capacity] : nullptr;
      capacity = other.capacity;
    }
    std::memcpy(bits, other.bits, capacity);
    return *this;
  }

  HuffmanBitCode &operator=(HuffmanBitCode &&other)
  {
    if (this == &other)
      return *this;

    capacity = other.capacity;
    length = other.length;
    byte = other.byte;
    std::swap(bits, other.bits);

    return *this;
  }

  HuffmanBitCode()
  {
  }

  ~HuffmanBitCode()
  {
    delete[] bits;
  }

  void binary_increment_bits()
  {
    bool carry = 1;
    for (int i = length - 1; i >= 0; i--)
    {
      bool bit = get_bit_at(bits[i / 8], i % 8);
      unsigned char result = carry + bit;
      set_bit_at(bits[i / 8], i % 8, result & 1U);
      carry = result & 2U;
    }
  }

  void append(bool bit)
  {
    if (capacity * 8 == length)
    {
      capacity++;
      unsigned char *new_bits = new unsigned char[capacity];
      std::memcpy(new_bits, bits, (length + 1) / 8);
      delete[] bits;
      bits = new_bits;
    }
    set_bit_at(bits[length / 8], length % 8, bit);
    length++;
  }

  void pop_bit()
  {
    length--;
  }
};

void write_byte(std::ostream &output, unsigned char byte)
{
  output.write((char *)&byte, sizeof(byte));
}

unsigned char read_byte(std::istream &input)
{
  unsigned char byte;
  input.read((char *)&byte, sizeof(byte));
  return byte;
}


void write_u32(std::ostream &output, uint32_t val)
{
  // TODO use ntohl
  output.write((char *)&val, sizeof(val));
}

uint32_t read_u32(std::istream &input)
{
  // TODO use ntohl
  uint32_t byte;
  input.read((char *)&byte, sizeof(byte));
  return byte;
}


void write_u64(std::ostream &output, uint64_t val)
{
  // TODO use ntohl
  output.write((char *)&val, sizeof(val));
}

uint64_t read_u64(std::istream &input)
{
  // TODO use ntohl
  uint64_t byte;
  input.read((char *)&byte, sizeof(byte));
  return byte;
}

struct HuffmanCodeBook
{
  std::vector<HuffmanBitCode> codes;
  std::vector<const HuffmanBitCode *> index = std::vector<const HuffmanBitCode *>(256);
  std::vector<unsigned char> length_increases_indices;
  int first_byte_idx;

  void fromTree(HuffmanNode *node, HuffmanBitCode &code)
  {
    if (node == nullptr)
    {
      // todo maybe not needed if, bin tree
      return;
    }

    if (node->children[0] == nullptr && node->children[1] == nullptr)
    {
      code.byte = node->byte;
      codes.push_back(code);
      return;
    }

    for (int i = 0; i < 2; i++)
    {
      code.append(i);
      fromTree(node->children[i], code);
      code.pop_bit();
    }
  }

  HuffmanCodeBook() {}

  HuffmanCodeBook(HuffmanNode *tree)
  {
    // construct non-canonical code from tree
    HuffmanBitCode base;
    fromTree(tree, base);

    // sort values first by length and then by alphabetical value
    std::sort(std::begin(codes), std::end(codes), [](const HuffmanBitCode &a, const HuffmanBitCode &b) {
      if (a.length != b.length)
      {
        return a.length < b.length;
      }

      for (int i = 0; i < a.length; i++)
      {
        bool a_bit = get_bit_at(a.bits[i / 8], i % 8);
        bool b_bit = get_bit_at(b.bits[i / 8], i % 8);
        if (a_bit != b_bit)
        {
          return a_bit < b_bit;
        }
      }
      return false;
    });

    // assign consecutive values using canonical huffman code algorithm
    HuffmanBitCode current_value;
    current_value.append(0);

    for (int i = 0; i < codes.size(); i++)
    {
      const HuffmanBitCode &code = codes[i];
      index[code.byte] = &code;

      if (current_value.length != code.length)
      {
        for (int j = current_value.length; j < code.length; j++)
        {
          length_increases_indices.push_back(i);
          current_value.append(0);
        }
      }

      std::memcpy(code.bits, current_value.bits, current_value.capacity); // TODO maybe change to length capacity, cause this is unsafe, if same length but diff cap

      current_value.binary_increment_bits();
    }
  }

  HuffmanNode *to_tree()
  {
    HuffmanNode *root = new HuffmanNode();
    for (const HuffmanBitCode &code : codes)
    {
      HuffmanNode *current = root;
      for (int j = 0; j < code.length; j++)
      {
        bool bit = get_bit_at(code.bits[j / 8], j % 8);

        if (current->children[bit] == nullptr)
        {
          current->children[bit] = new HuffmanNode();
        }

        if (j == code.length - 1)
        {
          current->children[bit]->byte = code.byte;
        }
        else
        {
          current = current->children[bit];
        }
      }
    }

    return root;
  }

  void write_to_file(std::ostream &output)
  {
    write_byte(output, length_increases_indices.size());
    for (int i = 0; i < length_increases_indices.size(); i++)
    {
      write_byte(output, length_increases_indices[i]);
    }

    // storing -1 because 0 size codebook is useless, and so 256 (max size) fits in 1 byte 
    write_byte(output, codes.size() - 1);
    for (int i = 0; i < codes.size(); i++)
    {
      write_byte(output, codes[i].byte);
    }
  }

  void read_from_file(std::istream &input)
  {
    int size = read_byte(input);
    for (int i = 0; i < size; i++)
    {
      int num = read_byte(input);

      length_increases_indices.push_back(num);
    }

    // assign consecutive values using canonical huffman code algorithm
    HuffmanBitCode current_value;
    current_value.append(0);

    size = read_byte(input) + 1;
    int length_inc_idx = 0;
    for (int i = 0; i < size; i++)
    {
      unsigned char byte = read_byte(input);

      while (length_inc_idx < length_increases_indices.size() && length_increases_indices[length_inc_idx] == i)
      {
        current_value.append(0);
        length_inc_idx++;
      }

      current_value.byte = byte;
      codes.push_back(current_value);
      index[byte] = &codes.back();

      current_value.binary_increment_bits();
    }
  }
};

void printCodeBook(const HuffmanCodeBook &b)
{
  for (const HuffmanBitCode &code : b.codes)
  {
    std::cout << (char)code.byte << ": ";
    for (int j = 0; j < code.length; j++)
    {
      std::cout << (int)get_bit_at(code.bits[j / 8], j % 8);
    }
    std::cout << std::endl;
  }
}

class HuffmanEncoder
{
public:
  HuffmanEncoder(std::istream &input, std::ostream &output) : input(input), output(output) {}
  void run()
  {
    build_tree();
    write_codebook();
    write_data();
  }

private:
  std::istream &input;
  std::ostream &output;

  HuffmanCodeBook codeBook;

  HuffmanNode *root;

  void build_tree()
  {
    std::vector<int> frequencies(256, 0);

    unsigned char byte;
    while (byte = read_byte(input), input)
    {
      frequencies[byte]++;
    }

    auto freq_cmp = [&](const HuffmanNode *a, const HuffmanNode *b) { return a->frequency > b->frequency; };

    std::priority_queue<HuffmanNode *, std::vector<HuffmanNode *>, decltype(freq_cmp)> queue(freq_cmp);

    for (unsigned i = 0; i < 256; i++)
    {
      unsigned char byte = i;
      if (frequencies[byte] > 0)
      {
        queue.push(new HuffmanNode(byte, frequencies[byte], nullptr, nullptr));
      }
    }

    // TODO check if queue is empty or maybe just one value

    while (queue.size() > 1)
    {
      // printBT(queue.top());
      // std::cout << "\n-----------\n";
      HuffmanNode *a = queue.top();
      queue.pop();
      // printBT(queue.top());
      // std::cout << "\n-----------\n";
      HuffmanNode *b = queue.top();
      queue.pop();
      queue.push(new HuffmanNode(0 /*whatever*/, a->frequency + b->frequency, a, b));
    }

    root = queue.top();
    //printBT(root);
    codeBook = HuffmanCodeBook(root);
    //printCodeBook(codeBook);
  }

  void write_codebook()
  {
    codeBook.write_to_file(output);
  }

  void write_data()
  {
    OBitStream encoded = OBitStream(output);

    input.clear();                 // clear fail and eof bits
    input.seekg(0, std::ios::beg); // back to the start!

    unsigned char byte;
    while (byte = read_byte(input), input)
    {
      const HuffmanBitCode *code = codeBook.index[byte];
      encoded.write_bit_chunk(code->bits, code->length);
    }
    encoded.flush();
  }
};

class HuffmanDecoder
{
public:
  HuffmanDecoder(std::istream &input, std::ostream &output) : input(input), output(output) {}
  void run()
  {
    read_codebook();
    build_tree();
    read_data();
  }

private:
  std::istream &input;
  std::ostream &output;

  HuffmanCodeBook codeBook;

  HuffmanNode *root;

  void read_codebook()
  {
    input.clear();                 // clear fail and eof bits
    input.seekg(0, std::ios::beg); // back to the start!

    codeBook.read_from_file(input);
    //printCodeBook(codeBook);
  }

  void build_tree()
  {
    root = codeBook.to_tree();

    printBT(root);
  }

  void read_data()
  {
    IBitStream reader(input);

    HuffmanNode *current = root;
    while (reader.has_more())
    {
      bool bit = reader.read_bit();
      if (current->children[bit] == nullptr)
      {
        write_byte(output, current->byte);
        current = root->children[bit];
      }
      else
      {
        current = current->children[bit];
      }
    }
    write_byte(output, current->byte);
  }
};

const uint64_t EMPTY_DIR_SPECIAL_VALUE = 0xFFFFFFFFFFFFFFFF;
struct FileEntry {
  uint64_t content_pointer;
  uint64_t content_length;
  uint32_t crc = 0xFFFFFFFF;
  uint64_t original_size;

  std::filesystem::path path;

  FileEntry() { }

  FileEntry(const std::filesystem::path& path) : path(path) {

  }

  void serialize(std::ostream& output) const {
    write_u64(output, content_pointer);
    write_u64(output, content_length);
    write_u32(output, crc);
    if (global_options & Option::has_compression_levels) {
      write_u64(output, original_size);
    }
    //strip root of path and normalize
    std::string filename = path.relative_path().lexically_normal().string();
    write_u32(output, (uint32_t)filename.size());

    for (char c : filename) {
      // TODO check if char to unsigned char conversion is safe
      write_byte(output, (unsigned char) c);
    }
  }

  void deserialize(std::istream& input) {
    content_pointer = read_u64(input);
    content_length = read_u64(input);
    crc = read_u32(input);

    if (global_options & Option::has_compression_levels) {
      original_size = read_u64(input);
    }

    std::string filename;
    uint32_t size = read_u32(input);
    for (int i = 0; i < size; i++) {
      // TODO check if unsigned char to char conversion is safe
      filename.push_back((char) read_byte(input));
    }
    path = filename;
  }
};

void assert_user(bool cond, std::string msg) {
    if (cond) {
        return;
    }
    // TODO throw
    std::cout << "ERROR: " << msg << std::endl;
    exit(1);
}

static uint32_t crc_table[256];

class Huffman {
  std::filesystem::path archive_path;
  std::fstream archive;
  std::fstream raw_output;
  std::vector<std::filesystem::path> file_paths;
  std::vector<FileEntry> fileEntries;
  uint64_t codebookPointer;
  uint32_t crc = 0xFFFFFFFF;
  std::vector<uint64_t> frequencies = std::vector<uint64_t>(256, 0);

  HuffmanCodeBook codeBook;
  HuffmanNode *tree;

public:
  Huffman(const std::filesystem::path& archive_path, const std::vector<std::filesystem::path>& file_paths, const std::filesystem::path& raw_output_path)
    : archive_path(archive_path),
      file_paths(file_paths),
      raw_output(raw_output_path, std::fstream::trunc | std::fstream::binary | std::fstream::in | std::fstream::out)
  {
    crc32::generate_table(crc_table);
  }

  void create()
  {
    fileEntries = std::vector<FileEntry>(file_paths.begin(), file_paths.end());

    archive.open(archive_path, std::fstream::trunc | std::fstream::binary | std::fstream::out | std::fstream::in);

    build_tree_from_files();  

    if (tree == nullptr) {
      global_options |= Option::has_only_empty_files;
    }

    write_header();
    write_files_content();
    write_codebook();
    write_frequencies();
    write_file_entries();

    write_u64(archive, codebookPointer);
    write_u32(archive, calculate_crc(archive.tellp()));
  }

  void extract()
  {
    archive.open(archive_path, std::fstream::binary | std::fstream::in);

    read_header();
    archive.seekg(-(64 + 32) / 8, std::ios::end);
    codebookPointer = read_u64(archive);
    crc = read_u32(archive);

    archive.seekg(codebookPointer, std::ios::beg);

    read_codebook();
    read_frequencies();
    build_tree();
    read_file_entries();
    extract_files();
  }

  void list()
  {
    archive.open(archive_path, std::fstream::binary | std::fstream::in);

    read_header();
    archive.seekg(-(64 + 32) / 8, std::ios::end);
    codebookPointer = read_u64(archive);
    crc = read_u32(archive);

    archive.seekg(codebookPointer, std::ios::beg);

    read_codebook();
    read_frequencies();
    read_file_entries();
    print_files();
  }

  void update()
  {
    archive.open(archive_path, std::fstream::binary | std::fstream::in | std::fstream::out);
    read_header();

    assert_user(global_options & Option::has_frequencies, "Create archive with -U to be able to update it fast");
    
    archive.seekg(-(64 + 32) / 8, std::ios::end);
    codebookPointer = read_u64(archive);
    crc = read_u32(archive);

    archive.seekg(codebookPointer, std::ios::beg);

    read_codebook();
    read_frequencies();
    
    read_file_entries();

    const FileEntry& update_file = file_paths[0];

    for (FileEntry& f : fileEntries) {
      if (f.path == update_file.path) {
        archive.seekp(codebookPointer, std::ios::beg);
        // add new byte frequencies to tree 
        build_tree_from_files({f});
        compress_file(f);
      }
    }
    
    write_codebook();
    write_frequencies();
    write_file_entries();

    write_u64(archive, codebookPointer);
    write_u32(archive, calculate_crc(archive.tellp()));
  }
  
  void check_integrity() {
    archive.open(archive_path, std::fstream::binary | std::fstream::in);
    read_header();
    archive.seekg(-(64 + 32) / 8, std::ios::end);
    codebookPointer = read_u64(archive);
    crc = read_u32(archive);
    archive.seekg(-4, std::ios::end);
    if (crc != calculate_crc(archive.tellg())) {
      std::cout << "File is corrupted" << std::endl;
    } else {
      std::cout << "File is probably not corrupted" << std::endl;
    }
  }

private:
  void write_header() {
    write_byte(archive, global_options);
  }

  void read_header() {
    global_options = read_byte(archive);
  }

  uint32_t calculate_crc(std::streampos to) {
    uint32_t crc = 0xFFFFFFFF;
    archive.clear();
    archive.seekg(0, std::ios::beg);

    unsigned char byte;
    while (archive.tellg() != to)
    {
      byte = read_byte(archive);
      crc = crc32::update(crc_table, crc, &byte, sizeof(byte));
      // std::cout << crc << std::endl;
    }
    return crc;
  }

  void read_frequencies() {
    if (global_options & Option::has_frequencies) {
      uint64_t size = read_u64(archive);
      for (int i = 0; i < size; i++)
      {
        frequencies[i] = read_u64(archive);
      }
    }
  }
  
  void write_frequencies() {
    if (global_options & Option::has_frequencies) {
      write_u64(archive, frequencies.size());
      for (int i = 0; i < frequencies.size(); i++)
      {
        write_u64(archive, frequencies[i]);
      }
    }
  }

  void build_tree_from_files() {
    build_tree_from_files(fileEntries);
  }

  void build_tree_from_files(const std::vector<FileEntry>& files)
  {
    unsigned char byte;
    for (const FileEntry& f : files) {
      std::ifstream input(f.path, std::fstream::binary); 
      while (byte = read_byte(input), input)
      {
        frequencies[byte]++;
      }
    }

    auto freq_cmp = [&](const HuffmanNode *a, const HuffmanNode *b) { return a->frequency > b->frequency; };

    std::priority_queue<HuffmanNode *, std::vector<HuffmanNode *>, decltype(freq_cmp)> queue(freq_cmp);

    for (unsigned i = 0; i < 256; i++)
    {
      unsigned char byte = i;
      if (frequencies[byte] > 0)
      {
        queue.push(new HuffmanNode(byte, frequencies[byte], nullptr, nullptr));
      }
    }

    // empty file
    if (queue.size() == 0) {
      queue.push(nullptr);
    } else if (queue.size() == 1) {
      HuffmanNode *a = queue.top();
      queue.pop();
      queue.push(new HuffmanNode(0, a->frequency, a, nullptr));
    } else {
      while (queue.size() > 1)
      {
        // printBT(queue.top());
        // std::cout << "\n-----------\n";
        HuffmanNode *a = queue.top();
        queue.pop();
        // printBT(queue.top());
        // std::cout << "\n-----------\n";
        HuffmanNode *b = queue.top();
        queue.pop();
        queue.push(new HuffmanNode(0 /*whatever*/, a->frequency + b->frequency, a, b));
      }
    }

    tree = queue.top();
    //printBT(tree);
    codeBook = HuffmanCodeBook(tree);
    //printCodeBook(codeBook);
  }

  void write_files_content()
  {
    for (FileEntry& f : fileEntries) {
      if (std::filesystem::is_directory(f.path)) {
        f.content_pointer = EMPTY_DIR_SPECIAL_VALUE;
      } else {
        compress_file(f);
      }
    }
  }

  void compress_file(FileEntry& f) {
    std::ifstream input(f.path, std::fstream::binary); 
    f.content_pointer = archive.tellg();

    OBitStream encoded = OBitStream(archive);

    int64_t input_size = 0;
    unsigned char byte;
    while (byte = read_byte(input), input)
    {
      const HuffmanBitCode *code = codeBook.index[byte];
      encoded.write_bit_chunk(code->bits, code->length);
      input_size++;
    }
    f.content_length = encoded.flush();

    f.original_size = input_size;
  }

  void write_codebook()
  {
    codebookPointer = archive.tellp();
    if (global_options & Option::has_only_empty_files) {
      return;
    }
    codeBook.write_to_file(archive);
  }

  void write_file_entries()
  { 
    write_u32(archive, fileEntries.size());

    for (const FileEntry& f : fileEntries) {
      f.serialize(archive);
    }
  }

  void read_codebook()
  {
    if (global_options & Option::has_only_empty_files) {
      return;
    }
    codeBook.read_from_file(archive);
    //std::cout << "Read codebook: " << std::endl;
    //printCodeBook(codeBook);
  }

  void build_tree()
  {
    tree = codeBook.to_tree();
    //std::cout << "Generated tree: " << std::endl;
    //printBT(tree);
  }
  
  void read_file_entries() {
    uint32_t size = read_u32(archive);

    for (int i = 0; i < size; i++) {
      FileEntry f;
      f.deserialize(archive);
      fileEntries.push_back(f);
    }
  }

  void extract_files() {
    for (FileEntry& f : fileEntries) {
      if (f.content_pointer == EMPTY_DIR_SPECIAL_VALUE) {
        std::filesystem::create_directories(f.path);
      } else {
        std::filesystem::path parent = f.path.parent_path();
        if (!parent.empty()) {
          std::filesystem::create_directories(parent);
        }
        std::ofstream output(f.path, std::fstream::trunc | std::fstream::binary);
        std::cout << "writing to " << f.path << std::endl;
        
        archive.seekp(f.content_pointer, std::ios::beg);
        IBitStream encoded = IBitStream(archive, f.content_length);

        HuffmanNode *current = tree;
        while (encoded.has_more())
        {
          bool bit = encoded.read_bit();
          if (current->children[bit] == nullptr)
          {
            write_byte(output, current->byte);
            current = tree->children[bit];
          }
          else
          {
            current = current->children[bit];
          }
        }
        write_byte(output, current->byte);
      }      
    }
  }

  void print_files() {
    std::cout.precision(2);
    for (FileEntry& f : fileEntries) {
        std::cout << f.path.string();

        if (global_options & Option::has_compression_levels) {
          std::cout << " ( "
            << (long long) f.original_size - (long long) f.content_length / 8 << " bytes saved; "
            << 100 - (float) f.content_length / 8 / f.original_size * 100 << "% )";
        }

        std::cout << std::endl;
    }
  }
};

/*
int main()
{
  std::stringstream stream;
  OBitStream writer(stream);

  writer.write_bit(1);
  writer.write_bit(0);
  writer.write_bit(1);
  writer.write_bit(0);

  writer.write_bit(1);
  writer.write_bit(0);
  writer.write_bit(1);
  writer.write_bit(0);

  uint64_t size = writer.flush();

  IBitStream reader(stream, size);

  std::cout << size << std::endl;

  while (reader.has_more()) {
    std::cout << (int)reader.read_bit() << std::endl;
  
  }

  
  // std::stringstream in("abccdddeeeeeffffffffggggggggggggghhhhhhhhhhhhhhhhhhhhhiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiijjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjjkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooopppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss");;
  // std::stringstream in("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Etiam ut suscipit erat. Integer massa urna, lacinia non ultricies non, eleifend eu enim. Etiam pretium turpis tellus, sed mollis metus dignissim at. Phasellus pellentesque magna quis convallis interdum. In a finibus velit, vel consequat ante. Curabitur blandit, eros quis laoreet accumsan, turpis turpis egestas felis, eget iaculis sapien felis ac mauris. Morbi ut nisl in ipsum tincidunt tincidunt lobortis eu libero. Quisque nibh risus, vehicula in justo at, sagittis commodo risus. Ut et neque in erat molestie aliquam. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Vivamus vulputate orci at molestie ultrices. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Maecenas pulvinar eleifend orci vel semper. Praesent ut mattis metus. Nulla sollicitudin nisi sit amet fermentum consequat. Donec malesuada tortor ut ligula lacinia ultricies.");
  std::stringstream in("1234567890qwertyuiopasdfghjklzxcvbnm");

  std::stringstream out;
  std::stringstream out2;

  HuffmanEncoder(in, out).run();
  
  IBitStream reader(out);
  while (reader.has_more()) {
    std::cout << (int)reader.read_bit();
  }
  std::cout << std::endl;
  std::cout << "----------" << std::endl;

  HuffmanDecoder(out, out2).run();

  std::cout << out2.str() << std::endl;

  std::cout << in.str().size() << " --> " << out.str().size() << " --> " << out2.str().size() << std::endl;
  

}
*/