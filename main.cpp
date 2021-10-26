#include <getopt.h>
#include <iostream>
#include <list>
#include <regex>

#include "huffman.h"

static struct option long_options[] = {
    // archive file
    {"file", required_argument, nullptr, 'f'},

    // add file, useful for files beginning with dash, or containing '*' and '?'
    {"add-file", required_argument, nullptr, 'a'},

    // change directory
    {"directory", required_argument, nullptr, 'C'},

    // raw mode output file (default stdout)
    {"output", required_argument, nullptr, 'o'},

    // only compress file content, not saving directory metadata
    {"raw", no_argument, nullptr, 'r'},

    // create archive from files
    {"create", no_argument, nullptr, 'c'},

    // extract files from archive
    {"extract", no_argument, nullptr, 'x'},

    // list archived files
    {"list", no_argument, nullptr, 'l'},

    // update file in archive
    {"update", no_argument, nullptr, 'u'},

    // check for archive integrity, and lists corrupted files
    {"check-integrity", no_argument, nullptr, 'i'},

    // adds integrity data to the archive
    {"add-integrity", no_argument, nullptr, 's'},

    // adds integrity data to the archive, and for each file
    {"add-integrity-each", no_argument, nullptr, 'e'},

    // adds copression levels to each file
    {"add-compression-levels", no_argument, nullptr, 'L'},

    // save frequencies to make updating faster
    {"save-frequencies", no_argument, nullptr, 'U'},

    {nullptr, 0, nullptr, 0}};

enum class Operation
{
    NONE,
    RAW,
    CREATE,
    EXTRACT,
    LIST,
    UPDATE,
    CHECK_INTEGRITY
};

std::regex glob_to_regex(const std::string& input) {
    // escape regex special characters except * and ?
    static std::regex regexSpecialChars { R"([-[\]{}()+.,\^$|#\s])" };
    static std::regex replaceQuestion { R"(\?)" };
    static std::regex replaceStar { R"(\*)" };

    std::string sanitized = std::regex_replace(input, regexSpecialChars, R"(\$&)" );

    std::string converted_q = std::regex_replace(sanitized, replaceQuestion, "." );
    std::string converted_s = std::regex_replace(converted_q, replaceStar, ".*" );

    // std::cout << "re: " << converted2 << std::endl;
    return std::regex(converted_s);
}

int main(int argc, char **argv)
{
    std::filesystem::path archive;
    std::filesystem::path change_dir;
    std::filesystem::path raw_output;
    std::vector<std::filesystem::path> files;
    Operation op = Operation::NONE;

    char ch;
    while ((ch = getopt_long(argc, argv, "f:a:C:o:rcxluiseLU", long_options, nullptr)) != -1)
    {
        switch (ch)
        {
        case 'f':
            assert_user(archive.empty(), "Specifying more than one archive file ('-f') not allowed");
            archive = optarg;
            break;
        case 'a':
            files.push_back(optarg);
            break;
        case 'C':
            assert_user(change_dir.empty(), "Multiple -C not supported.");
            change_dir = optarg;
            break;
        case 'o':
            raw_output = optarg;
            break;
        case 'r':
            assert_user(op == Operation::NONE, "Specifying more than one of '-rcxlui' not allowed");
            op = Operation::RAW;
            break;
        case 'c':
            assert_user(op == Operation::NONE, "Specifying more than one of '-rcxlui' not allowed");
            op = Operation::CREATE;
            break;
        case 'x':
            assert_user(op == Operation::NONE, "Specifying more than one of '-rcxlui' not allowed");
            op = Operation::EXTRACT;
            break;
        case 'l':
            assert_user(op == Operation::NONE, "Specifying more than one of '-rcxlui' not allowed");
            op = Operation::LIST;
            break;
        case 'u':
            assert_user(op == Operation::NONE, "Specifying more than one of '-rcxlui' not allowed");
            op = Operation::UPDATE;
            break;
        case 'i':
            assert_user(op == Operation::NONE, "Specifying more than one of '-rcxlui' not allowed");
            op = Operation::CHECK_INTEGRITY;
            break;
        case 's':
            global_options |= Option::has_crc;
            break;
        case 'e':
            global_options |= Option::has_crc_each;
            break;
        case 'L':
            global_options |= Option::has_compression_levels;
            break;
        case 'U':
            global_options |= Option::has_frequencies;
            break;
        default:
            std::cout << "Usage: " << std::endl;
            return EXIT_FAILURE;
        }
    }

    assert_user(!archive.empty(), "You must specify an archive file name ('-f')");

    archive = std::filesystem::absolute(archive);

    if (!raw_output.empty())
    {
        raw_output = std::filesystem::absolute(archive);
    }

    if (!change_dir.empty())
    {
        assert_user(std::filesystem::exists(change_dir), "Change dir " + change_dir.string() + " does not exist");
        std::filesystem::current_path(change_dir);
    }

    for (int i = optind; i < argc; i++)
    {
        std::string filename = argv[i];
        std::filesystem::path path = filename;

        std::vector<std::filesystem::path> reached;
        if (filename.find('*') != std::string::npos || filename.find('?') != std::string::npos)
        {
            reached = {""};
            std::vector<std::filesystem::path> new_reached;
            
            for (const std::filesystem::path& part : path) {
                std::regex re = glob_to_regex(part);
                std::string path_string = part.string();
                // TODO move this up
                for (const std::filesystem::path& curr : reached) {
                    if (
                        path_string.find('*') != std::string::npos ||
                        path_string.find('?') != std::string::npos
                    )
                    {
                        if (!std::filesystem::is_directory(curr)) {
                            continue;
                        }
                        for (const std::filesystem::directory_entry& e: std::filesystem::directory_iterator(curr)) {
                            std::string last_filename = e.path().filename();


                            if (std::regex_match(last_filename, re)) {
                                new_reached.push_back(curr / last_filename);
                            }
                        }
                    }
                    else
                    {
                        new_reached.push_back(curr / part);
                    }
                }
                // for (const std::filesystem::path& r : new_reached) {
                //     std::cout << r << std::endl;
                // }
                // std::cout << "-----------------" << std::endl;
                reached = std::move(new_reached);
                new_reached = {};
            }
            for (const std::filesystem::path& r : reached) {
                if (std::filesystem::exists(r)) {
                    new_reached.push_back(r);
                }
            }
            reached = std::move(new_reached);
        }
        else 
        {
            reached = {path};
        }

        for (const std::filesystem::path& reached_path : reached) {
            if (std::filesystem::is_directory(reached_path))
            {
                if (std::filesystem::is_empty(reached_path))
                {
                    files.push_back(reached_path);
                }
                else
                {
                    for (const std::filesystem::directory_entry& e : std::filesystem::recursive_directory_iterator(reached_path))
                    {
                        files.push_back(e.path());
                    }
                }
            }
            else
            {
                files.push_back(reached_path);
            }
        }
    }

    // TODO set files

    for (const std::filesystem::path &f : files)
    {
        assert_user(std::filesystem::exists(f), "File " + f.string() + " does not exist");
    }

    assert_user(op != Operation::NONE, "You must specify one of '-rcxlui' operation");

    Huffman huffman(archive, files, raw_output);

    switch (op)
    {
    case Operation::CREATE:
        assert_user(!files.empty(), "Refusing to create an empty archive");
        huffman.create();
        break;
    case Operation::EXTRACT:
        assert_user(std::filesystem::exists(archive), "File " + archive.string() + " does not exist");
        huffman.extract();
        break;
    case Operation::LIST:
        assert_user(std::filesystem::exists(archive), "File " + archive.string() + " does not exist");
        huffman.list();
        break;
    case Operation::UPDATE:
        assert_user(files.size() == 1, "Please specify a file to update");
        huffman.update();
        break;
    case Operation::CHECK_INTEGRITY:
        assert_user(std::filesystem::exists(archive), "File " + archive.string() + " does not exist");
        huffman.check_integrity();
        break;
    }
}