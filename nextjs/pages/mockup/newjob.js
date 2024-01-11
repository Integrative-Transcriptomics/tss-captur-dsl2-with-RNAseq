import React from 'react';
import Layout from '../components/mockup/Layout';

export default function FormPage() {
    const handleSubmit = async (event) => {
        // submit
    };

    return (
        <Layout>
            <form onSubmit={handleSubmit}>
                <label>
                    TSV File:
                    <input type="file" name="tsvFile" accept=".tsv" />
                </label>
                <label>
                    Input Genomes:
                    <input type="file" name="inputGenomes" webkitdirectory="" directory="" />
                </label>
                <label>
                    Input GFFs:
                    <input type="file" name="inputGFFs" webkitdirectory="" directory="" />
                </label>
                <input type="submit" value="Submit" />
            </form>
        </Layout>
    );
}