#pragma once
#include <string>
#include <string_view>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <iostream> // Added for logging bytes
#include <iomanip>  // Added for hex formatting

struct PublicKeyEntry {
    std::string groupId;
    std::vector<unsigned char> pubKeyBytes;
};

class Helper {
public:
    // 1. CD-Key Base-24 Encoder (Extract CD key bytes) with direct screen print for each byte
    static std::vector<unsigned char> EncodeBinaryKey(const std::string& cdKey) {
        std::string alphabet = "BCDFGHJKMPQRTVWXY2346789";
        std::string rawKey = "";
        for (char c : cdKey) {
            if (c != '-') rawKey += (char)std::toupper(c);
        }

        if (rawKey.length() != 25) {
            throw std::runtime_error("Key must be 25 characters.");
        }

        std::vector<unsigned char> digits(25, 0);
        bool isNKey = false;
        int digitCount = 0;

        for (char ch : rawKey) {
            if (ch == 'N' && !isNKey) {
                isNKey = true;
                for (int i = digitCount; i > 0; i--) {
                    digits[i] = digits[i - 1];
                }
                digits[0] = (unsigned char)digitCount;
                digitCount++;
                continue;
            }

            size_t val = alphabet.find(ch);
            if (val == std::string::npos) {
                throw std::runtime_error("Invalid character in key.");
            }
            digits[digitCount] = (unsigned char)val;
            digitCount++;
        }

        std::vector<unsigned char> binary(16, 0);
        for (int digit : digits) {
            uint32_t carry = digit;
            for (int i = 0; i < 16; i++) {
                uint32_t res = (binary[i] * 24) + carry;
                binary[i] = (unsigned char)(res & 0xFF);
                carry = res >> 8;
            }
        }

        if (isNKey) {
            binary[14] |= 0x08;
        }

        return binary;
    }

    // 2. Extract Group and PublicKey entries from XML/Base64 configuration data with byte logging
    static bool ParseKeyEntriesSimple(const std::string& outerXml, std::vector<PublicKeyEntry>& pkEntries,
        std::vector<unsigned char>(*base64DecodeFunc)(std::string_view)) {
        size_t infoBinPos = outerXml.find("pkeyConfigData");
        if (infoBinPos == std::string::npos) {
            std::cout << "[Debug] 'pkeyConfigData' tag not found in XML.\n";
            return false;
        }

        size_t contentStart = outerXml.find('>', infoBinPos) + 1;
        size_t contentEnd = outerXml.find("</", contentStart);
        std::string_view base64InfoBin(&outerXml[contentStart], contentEnd - contentStart);

        std::string cleanedB64;
        cleanedB64.reserve(base64InfoBin.size());
        for (char c : base64InfoBin) {
            if (!std::isspace((unsigned char)c)) cleanedB64.push_back(c);
        }

        std::vector<unsigned char> decodedInnerBytes = base64DecodeFunc(cleanedB64);
        if (decodedInnerBytes.empty()) {
            std::cout << "[Debug] Failed to decode inner Base64 block.\n";
            return false;
        }

        std::string innerXml(decodedInnerBytes.begin(), decodedInnerBytes.end());
        if (decodedInnerBytes.size() > 2 && decodedInnerBytes[1] == 0) {
            const wchar_t* wstr_ptr = reinterpret_cast<const wchar_t*>(decodedInnerBytes.data());
            int wstr_len = static_cast<int>(decodedInnerBytes.size() / 2);

            int sz = WideCharToMultiByte(CP_UTF8, 0, wstr_ptr, wstr_len, NULL, 0, NULL, NULL);
            std::string utf8Str(sz, 0);
            WideCharToMultiByte(CP_UTF8, 0, wstr_ptr, wstr_len, &utf8Str[0], sz, NULL, NULL);
            innerXml = utf8Str;
        }

        std::string_view innerSv = innerXml;
        size_t searchPos = 0;
        int blockIndex = 0;

        while (true) {
            size_t pkStart = innerSv.find("<pkc:PublicKey>", searchPos);
            if (pkStart == std::string_view::npos) {
                pkStart = innerSv.find("<PublicKey>", searchPos);
                if (pkStart == std::string_view::npos) break;
            }

            size_t pkEnd = innerSv.find("</pkc:PublicKey>", pkStart);
            if (pkEnd == std::string_view::npos) {
                pkEnd = innerSv.find("</PublicKey>", pkStart);
                if (pkEnd == std::string_view::npos) break;
                pkEnd += sizeof("</PublicKey>") - 1;
            }
            else {
                pkEnd += sizeof("</pkc:PublicKey>") - 1;
            }

            std::string_view block = innerSv.substr(pkStart, pkEnd - pkStart);
            searchPos = pkEnd;

            std::string groupId = "";
            size_t gStart = block.find("GroupId>");
            if (gStart == std::string_view::npos) gStart = block.find(":GroupId>");
            if (gStart != std::string_view::npos) {
                gStart = block.find('>', gStart) + 1;
                size_t gEnd = block.find("</", gStart);
                groupId = std::string(block.substr(gStart, gEnd - gStart));
            }

            size_t kStart = block.find("PublicKeyValue>");
            if (kStart == std::string_view::npos) kStart = block.find(":PublicKeyValue>");
            if (kStart != std::string_view::npos) {
                kStart = block.find('>', kStart) + 1;
                size_t kEnd = block.find("</", kStart);
                std::string_view rawKeyB64 = block.substr(kStart, kEnd - kStart);

                std::string cleanKeyB64;
                cleanKeyB64.reserve(rawKeyB64.size());
                for (char c : rawKeyB64) {
                    if (!std::isspace((unsigned char)c)) cleanKeyB64.push_back(c);
                }

                std::vector<unsigned char> pubKeyBytes = base64DecodeFunc(cleanKeyB64);

                if (pubKeyBytes.size() == 1579) {
                    pkEntries.push_back({ groupId, pubKeyBytes });
                }
            }
        }

        return !pkEntries.empty();
    }
};