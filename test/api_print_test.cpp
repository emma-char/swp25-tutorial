/ SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>

#include <print.hpp>

#include "app_test.hpp"

// To prevent issues when running multiple API tests in parallel, give each API test unique names:
struct print_test : public app_test
{};

TEST_F(print_test, out_empty)
{
    ::testing::internal::CaptureStdout();
    print(2);
    print("Hello!");
    std::string const std_out = testing::internal::GetCapturedStdout();

    std::string const expected = "2\nHello\n";

    EXPECT_EQ(std_out, expected);
    
}