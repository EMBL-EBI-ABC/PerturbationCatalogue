import { Component, OnInit } from '@angular/core';
import { HeaderComponent } from "../shared/header/header.component";
import { MatCardModule } from "@angular/material/card";
import { FooterComponent } from "../shared/footer/footer.component";
import { ElasticService } from '../services/elastic.service';
import { MatTableModule } from '@angular/material/table';
import { MatPaginator, MatPaginatorModule, PageEvent } from '@angular/material/paginator';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MatFormFieldModule } from '@angular/material/form-field';
import { MatInputModule } from '@angular/material/input';
import { marked } from 'marked';

@Component({
  selector: 'app-home',
  standalone: true,
  imports: [
    HeaderComponent,
    FooterComponent,
    MatCardModule,
    MatTableModule,
    MatPaginator,
    MatPaginatorModule,
    CommonModule,
    FormsModule,
    MatFormFieldModule,
    MatInputModule
  ],
  templateUrl: './home.component.html',
  styleUrl: './home.component.css'
})
export class HomeComponent implements OnInit {
  data: any;
  totalResults = 0;
  pageSize = 15;
  currentPage = 0;
  searchQuery: string = '';

  constructor(private elasticService: ElasticService) {}

  ngOnInit() {
    this.fetchData();
  }

  getKeys(obj: any): string[] {
    return Object.keys(obj);
  }

  fetchData() {
    const start = this.currentPage * this.pageSize;
    const size = this.pageSize;
    const query = this.searchQuery;
    this.elasticService.getData(start, size, query).subscribe(
      (response) => {
        this.data = response.results.map((row: Record<string, string>) => {
          const parsedRow: Record<string, string> = {};
          for (const key of this.getKeys(row)) {
            parsedRow[key] = marked.parse(row[key].toString()) as string;
          }
          return parsedRow;
        });
        this.totalResults = response.total;
      },
      (error) => {
        console.error('Error fetching data:', error);
      }
    );
  }

  onPageChange(event: PageEvent) {
    this.pageSize = event.pageSize;
    this.currentPage = event.pageIndex;
    this.fetchData();
  }

  onSearchChange(): void {
    this.currentPage = 0;
    this.fetchData();
  }
}
